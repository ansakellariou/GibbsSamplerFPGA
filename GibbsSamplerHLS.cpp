#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#define T 10
#define K 15
#define MOTIF_LENGTH 1973

unsigned int myseed = 0x015A4E36;

void profile_with_pseudocounts_without_i(char *motifs, int pseudocount, int motif_length, int number_of_motifs, int i_to_exclude, int *profile_matrix)
{
  #pragma HLS INLINE
  int i, j, l;
  int count_a, count_c, count_g, count_t, sum;
  // ******************************* profile_matrix ********************************** //
  // matrix of size 4 * motif_length ; each line corresponds to a single nucleotide    //
  // in the order A, C, G, T ; each column corresponds to the number of appearances    //
  // of a certain nucleotide (depends on the line) in the same column in motifs matrix //
  // in this number we add the pseudocount to prevent zeros in the matrix and divide   //
  // the result with 2 * number_of_motifs to obtain the corresponding probability      //

  // initialization of profile_matrix //

  loop_profile_with_pseudocounts_init: for (i = 0 ; i < 4 ; i++)
  {
	#pragma HLS loop_tripcount min=4 max=4
    for (j = 0; j < K; j++)
    {
      #pragma HLS loop_tripcount min=15 max=15
      *((profile_matrix + i*4) + j) = pseudocount;
    }
  }


  loop_profile_with_pseudocounts: for (j = 0; j < K; j++)
  {
    #pragma HLS loop_tripcount min=15 max=15
    count_a = 0;
    count_c = 0;
    count_g = 0;
    count_t = 0;

    for(i = 0; i < T; i++)
    {
      #pragma HLS loop_tripcount min=10 max=10

      if (i == i_to_exclude)
        continue;
      if((*((motifs+i*number_of_motifs)+j)) == 'A')
        count_a++;
      else if((*((motifs+i*number_of_motifs)+j)) == 'C')
        count_c++;
      else if((*((motifs+i*number_of_motifs)+j)) == 'G')
        count_g++;
      else if((*((motifs+i*number_of_motifs)+j)) == 'T')
        count_t++;
    }

    sum = count_a + count_c + count_g + count_t + 4*pseudocount;

    *((profile_matrix + 0*4) + j) = *((profile_matrix + 0*4) + j) + count_a;
    *((profile_matrix + 1*4) + j) = *((profile_matrix + 1*4) + j) + count_c;
    *((profile_matrix + 2*4) + j) = *((profile_matrix + 2*4) + j) + count_g;
    *((profile_matrix + 3*4) + j) = *((profile_matrix + 3*4) + j) + count_t;
  }

  return;
}

int score(char *motifs, int motif_length, int number_of_motifs)
{
  #pragma HLS INLINE
  int i, j;
  int count_a, count_c, count_g, count_t;
  int score; // the number of unpopular (lower case) letters in the motif matrix //

  score = 0;

  // parallelizable //
  loop_score: for (j = 0; j < K; j++)
  {
    #pragma HLS loop_tripcount min=15 max=15
    count_a = 0;
    count_c = 0;
    count_g = 0;
    count_t = 0;

    for(i = 0; i < T; i++)
    {
	  #pragma HLS loop_tripcount min=10 max=10
      if((*((motifs+i*number_of_motifs)+j)) == 'A')
        count_a++;
      else if((*((motifs+i*number_of_motifs)+j)) == 'C')
        count_c++;
      else if((*((motifs+i*number_of_motifs)+j)) == 'G')
        count_g++;
      else if((*((motifs+i*number_of_motifs)+j)) == 'T')
        count_t++;
    }

    if ((count_a >= count_c) && (count_a >= count_g) && (count_a >= count_t))
      score = score + count_c + count_g + count_t;
    else if((count_c >= count_a) && (count_c >= count_g) && (count_c >= count_t))
      score = score + count_a + count_g + count_t;
    else if((count_g >= count_a) && (count_g >= count_c) && (count_g >= count_t))
      score = score + count_a + count_c + count_t;
    else if((count_t >= count_a) && (count_t >= count_c) && (count_t >= count_g))
      score = score + count_a + count_c + count_g;
  }

  return score;
}
void generate_probabilities_of_all_kmers(char *motif_excluded, int k, int motif_length, int t, int *profile_matrix, int *probabilities_matrix)
{
  #pragma HLS INLINE
  int i, j, temp;

  /*loop_generate_probabilities_of_all_kmers_init: for (i = 0; i < MOTIF_LENGTH - k + 1 ; i++)
  {
    #pragma HLS loop_tripcount min=1959 max=1959
    temp = 1;
  }*/
  loop_generate_probabilities_of_all_kmers: for (i = 0; i < MOTIF_LENGTH - K + 1 ; i++)
  {
    #pragma HLS loop_tripcount min=1959 max=1959
	  *(probabilities_matrix + i) = 1;
    // parallelizable //
    for (j = 0; j < K; j++)
    {
	  #pragma HLS loop_tripcount min=15 max=15
      if (*(motif_excluded+i+j) == 'A')
    	  *(probabilities_matrix + i)  = *(probabilities_matrix + i)  * (*((profile_matrix ) + j)) ; //0*4
      else if (*(motif_excluded+i+j) == 'C')
    	  *(probabilities_matrix + i)  = *(probabilities_matrix + i)  * (*((profile_matrix + 4) + j)); //1*4
      else if (*(motif_excluded+i+j) == 'G')
    	  *(probabilities_matrix + i)  = *(probabilities_matrix + i)  * (*((profile_matrix + 8) + j)); //2*4
      else if (*(motif_excluded+i+j) == 'T')
    	  *(probabilities_matrix + i)  = *(probabilities_matrix + i)  * (*((profile_matrix + 12) + j)); //3*4
    }
    //*(probabilities_matrix + i) = temp;
  }
  return;
}

void generate_random_kmer(char *motif_excluded, int k, int motif_length, int t, int *profile_matrix, char *most_probable_kmer)
{
    #pragma HLS INLINE
	int probabilities_matrix[MOTIF_LENGTH - K + 1];
    #pragma HLS ARRAY_PARTITION variable=probabilities_matrix block factor=4 dim=1
	int most_probable_i, i;
	int most_probable_val;
	int unfair_dice_result;

	generate_probabilities_of_all_kmers(motif_excluded, K, MOTIF_LENGTH, T, (int*)profile_matrix, probabilities_matrix);

	most_probable_i = 0;
	most_probable_val = 0;

	loop_generate_random_kmer: for (i = 0; i < MOTIF_LENGTH - K + 1 ; i++)
	{
		#pragma HLS loop_tripcount min=1959 max=1959
		if ((*(probabilities_matrix + i)) >= most_probable_val)
		{
		  most_probable_val = (*(probabilities_matrix + i));
		  most_probable_i = i;
		}
	}
	most_probable_kmer[0] = *(motif_excluded + most_probable_i+0);
	most_probable_kmer[1] = *(motif_excluded + most_probable_i+1);
	most_probable_kmer[2] = *(motif_excluded + most_probable_i+2);
	most_probable_kmer[3] = *(motif_excluded + most_probable_i+3);
	most_probable_kmer[4] = *(motif_excluded + most_probable_i+4);
	most_probable_kmer[5] = *(motif_excluded + most_probable_i+5);
	most_probable_kmer[6] = *(motif_excluded + most_probable_i+6);
	most_probable_kmer[7] = *(motif_excluded + most_probable_i+7);
	most_probable_kmer[8] = *(motif_excluded + most_probable_i+8);
	most_probable_kmer[9] = *(motif_excluded + most_probable_i+9);
	most_probable_kmer[10] = *(motif_excluded + most_probable_i+10);
	most_probable_kmer[11] = *(motif_excluded + most_probable_i+11);
	most_probable_kmer[12] = *(motif_excluded + most_probable_i+12);
	most_probable_kmer[13] = *(motif_excluded + most_probable_i+13);
	most_probable_kmer[14] = *(motif_excluded + most_probable_i+14);
	most_probable_kmer[15] = '\0';
  /*string_n_copy(most_probable_kmer, (motif_excluded + most_probable_i), k*sizeof(char));
  most_probable_kmer[k] = '\0';*/
  return;
}

int myrand()
{
  #pragma HLS INLINE
  unsigned int t = myseed * 0x015A4E35 + 1;
  myseed = t;
  return (int)(t >> 16) & 0x7FFF;
}


extern "C"
{
  /**
   * Extern is a requirement for Vitis Unified Software Platform,
   * when we use a .cpp (C++ source code) file.
   * */

  /***********************************************************
  * Function:  gibbsSamplerKernel
  *************************************************************/

  void gibbsSamplerKernel(char dna_input[T*(MOTIF_LENGTH+1)], char best_motifs_output[T*(K+1)], int k, int t, int n, int motif_length)
  {
    #pragma HLS INTERFACE s_axilite port=return bundle=control

    #pragma HLS INTERFACE m_axi port=best_motifs_output offset=slave bundle=gmem
    #pragma HLS INTERFACE s_axilite port=best_motifs_output    bundle=control

    #pragma HLS INTERFACE m_axi port=dna_input offset=slave bundle=gmem
    #pragma HLS INTERFACE s_axilite port=dna_input    bundle=control

    /*** Required INTERFACE pragma END ***/

    #pragma HLS INTERFACE s_axilite port=k bundle=control
    #pragma HLS INTERFACE s_axilite port=t bundle=control
    #pragma HLS INTERFACE s_axilite port=n bundle=control
    #pragma HLS INTERFACE s_axilite port=motif_length bundle=control

    #pragma HLS INLINE
    int i, j, l, m, starting_position_of_kmer, i_to_exclude;
    char dna[T*(MOTIF_LENGTH+1)];
    #pragma HLS ARRAY_PARTITION variable=dna block factor=4 dim=1
    int profile_matrix[4][K];
    #pragma HLS ARRAY_PARTITION variable=profile_matrix block factor=4 dim=1
    char most_probable_kmer[K+1];
    #pragma HLS ARRAY_PARTITION variable=most_probable_kmer complete dim=1
    int score_of_motifs, score_of_best_motifs;
    char motifs[T][K+1];
    #pragma HLS ARRAY_PARTITION variable=motifs complete dim=1
    char best_motifs[T][K+1];
    #pragma HLS ARRAY_PARTITION variable=best_motifs complete dim=1
    int random_val, i_mul_k;

    memcpy(dna, dna_input,sizeof(char)*(T*(MOTIF_LENGTH+1)));

    printf("Motifs: \n");
    loop_gibbs_sampler_init: for (i = 0 ; i < T ; i++)
    {
    #pragma HLS loop_tripcount min=10 max=10
      random_val = myrand();
      starting_position_of_kmer = abs(random_val % (MOTIF_LENGTH - K + 1));

    motifs[i][0] = *(dna + i*(MOTIF_LENGTH+1) + starting_position_of_kmer+0);
    motifs[i][1] = *(dna + i*(MOTIF_LENGTH+1) + starting_position_of_kmer+1);
    motifs[i][2] = *(dna + i*(MOTIF_LENGTH+1) + starting_position_of_kmer+2);
    motifs[i][3] = *(dna + i*(MOTIF_LENGTH+1) + starting_position_of_kmer+3);
    motifs[i][4] = *(dna + i*(MOTIF_LENGTH+1) + starting_position_of_kmer+4);
    motifs[i][5] = *(dna + i*(MOTIF_LENGTH+1) + starting_position_of_kmer+5);
    motifs[i][6] = *(dna + i*(MOTIF_LENGTH+1) + starting_position_of_kmer+6);
    motifs[i][7] = *(dna + i*(MOTIF_LENGTH+1) + starting_position_of_kmer+7);
    motifs[i][8] = *(dna + i*(MOTIF_LENGTH+1) + starting_position_of_kmer+8);
    motifs[i][9] = *(dna + i*(MOTIF_LENGTH+1) + starting_position_of_kmer+9);
    motifs[i][10] = *(dna + i*(MOTIF_LENGTH+1) + starting_position_of_kmer+10);
    motifs[i][11] = *(dna + i*(MOTIF_LENGTH+1) + starting_position_of_kmer+11);
    motifs[i][12] = *(dna + i*(MOTIF_LENGTH+1) + starting_position_of_kmer+12);
    motifs[i][13] = *(dna + i*(MOTIF_LENGTH+1) + starting_position_of_kmer+13);
    motifs[i][14] = *(dna + i*(motif_length+1) + starting_position_of_kmer+14);
    motifs[i][15] = '\0';
      /*string_n_copy(motifs[i], (dna + i*(motif_length+1) + starting_position_of_kmer), k);
      motifs[i][k] = '\0';*/

    best_motifs[i][0] = motifs[i][0];
    best_motifs[i][1] = motifs[i][1];
    best_motifs[i][2] = motifs[i][2];
    best_motifs[i][3] = motifs[i][3];
    best_motifs[i][4] = motifs[i][4];
    best_motifs[i][5] = motifs[i][5];
    best_motifs[i][6] = motifs[i][6];
    best_motifs[i][7] = motifs[i][7];
    best_motifs[i][8] = motifs[i][8];
    best_motifs[i][9] = motifs[i][9];
    best_motifs[i][10] = motifs[i][10];
    best_motifs[i][11] = motifs[i][11];
    best_motifs[i][12] = motifs[i][12];
    best_motifs[i][13] = motifs[i][13];
    best_motifs[i][14] =  motifs[i][14];
    best_motifs[i][15] = '\0';
      /*string_n_copy(best_motifs[i], motifs[i], k);
      best_motifs[i][k] = '\0';*/
      printf("%s\n", motifs[i]);
    }


    /*for (i = 0 ; i < t ; i++)
    {
      #pragma HLS loop_tripcount min=10 max=10
      printf("%s\n", motifs[i]);
    }*/

    loop_gibbs_sampler: for (j = 0; j < n; j++)
    {
    #pragma HLS loop_tripcount min=500 max=500
      random_val = myrand();
      i_to_exclude = abs(random_val % t);

      profile_with_pseudocounts_without_i((char*)motifs, 1, K, T, i_to_exclude, (int*)profile_matrix);

      generate_random_kmer(dna + i_to_exclude*(MOTIF_LENGTH+1), K, MOTIF_LENGTH, t, (int*)profile_matrix, most_probable_kmer);

      *(*(motifs + i_to_exclude) + 0) = most_probable_kmer[0] ;
    *(*(motifs + i_to_exclude) + 1) = most_probable_kmer[1];
    *(*(motifs + i_to_exclude) + 2) = most_probable_kmer[2];
    *(*(motifs + i_to_exclude) + 3) = most_probable_kmer[3];
    *(*(motifs + i_to_exclude) + 4) = most_probable_kmer[4];
    *(*(motifs + i_to_exclude) + 5) = most_probable_kmer[5];
    *(*(motifs + i_to_exclude) + 6) = most_probable_kmer[6];
    *(*(motifs + i_to_exclude) + 7) = most_probable_kmer[7];
    *(*(motifs + i_to_exclude) + 8) = most_probable_kmer[8];
    *(*(motifs + i_to_exclude) + 9) = most_probable_kmer[9];
    *(*(motifs + i_to_exclude) + 10) = most_probable_kmer[10];
    *(*(motifs + i_to_exclude) + 11) = most_probable_kmer[11];
    *(*(motifs + i_to_exclude) + 12) = most_probable_kmer[12];
    *(*(motifs + i_to_exclude) + 13) = most_probable_kmer[13];
    *(*(motifs + i_to_exclude) + 14) = most_probable_kmer[14];
    *(*(motifs + i_to_exclude) + 15) = '\0';
      /*string_n_copy(*(motifs + i_to_exclude), most_probable_kmer, k);
      *(*(motifs + i_to_exclude) + k) = '\0';*/

      score_of_motifs = score((char*)motifs, K, T);
      score_of_best_motifs = score((char*)best_motifs, K, T);

      if(score_of_motifs < score_of_best_motifs)
      {
        *(*(best_motifs + i_to_exclude) + 0) = most_probable_kmer[0] ;
      *(*(best_motifs + i_to_exclude) + 1) = most_probable_kmer[1];
      *(*(best_motifs + i_to_exclude) + 2) = most_probable_kmer[2];
      *(*(best_motifs + i_to_exclude) + 3) = most_probable_kmer[3];
      *(*(best_motifs + i_to_exclude) + 4) = most_probable_kmer[4];
      *(*(best_motifs + i_to_exclude) + 5) = most_probable_kmer[5];
      *(*(best_motifs + i_to_exclude) + 6) = most_probable_kmer[6];
      *(*(best_motifs + i_to_exclude) + 7) = most_probable_kmer[7];
      *(*(best_motifs + i_to_exclude) + 8) = most_probable_kmer[8];
      *(*(best_motifs + i_to_exclude) + 9) = most_probable_kmer[9];
      *(*(best_motifs + i_to_exclude) + 10) = most_probable_kmer[10];
      *(*(best_motifs + i_to_exclude) + 11) = most_probable_kmer[11];
      *(*(best_motifs + i_to_exclude) + 12) = most_probable_kmer[12];
      *(*(best_motifs + i_to_exclude) + 13) = most_probable_kmer[13];
      *(*(best_motifs + i_to_exclude) + 14) = most_probable_kmer[14];
      *(*(best_motifs + i_to_exclude) + 15) = '\0';
        /*string_n_copy(*(best_motifs + i_to_exclude), most_probable_kmer, k);
        *(*(best_motifs + i_to_exclude) + k) = '\0';*/
      }
    }
    memcpy(best_motifs_output, best_motifs,sizeof(char)*T*(K+1));

    /*for(i = 0 ; i < T ; i++)
    {
    #pragma HLS loop_tripcount min=10 max=10
      i_mul_k = i*(K+1);
      *(best_motifs_output + i_mul_k + 0) = best_motifs[i][0] ;
    *(best_motifs_output + i_mul_k + 1) = best_motifs[i][1];
    *(best_motifs_output + i_mul_k + 2) = best_motifs[i][2];
    *(best_motifs_output + i_mul_k + 3) = best_motifs[i][3];
    *(best_motifs_output + i_mul_k + 4) = best_motifs[i][4];
    *(best_motifs_output + i_mul_k + 5) = best_motifs[i][5];
    *(best_motifs_output + i_mul_k + 6) = best_motifs[i][6];
    *(best_motifs_output + i_mul_k + 7) = best_motifs[i][7];
    *(best_motifs_output + i_mul_k + 8) = best_motifs[i][8];
    *(best_motifs_output + i_mul_k + 9) = best_motifs[i][9];
    *(best_motifs_output + i_mul_k + 10) = best_motifs[i][10];
    *(best_motifs_output + i_mul_k + 11) = best_motifs[i][11];
    *(best_motifs_output + i_mul_k + 12) = best_motifs[i][12];
    *(best_motifs_output + i_mul_k + 13) = best_motifs[i][13];
    *(best_motifs_output + i_mul_k + 14) = best_motifs[i][14];
    *(best_motifs_output + i_mul_k + 15) = '\0';
      /*string_n_copy(best_motifs_output+i*(k+1), *(best_motifs + i), k);
      *(best_motifs_output + i*(k+1) + k) = '\0';*/
    //}*/

    //score_of_motifs = score((char*)best_motifs, K, T);
    printf("Score of best motifss: %d\n", score_of_best_motifs);
    return;
  }

}
