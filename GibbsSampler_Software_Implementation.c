#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

char **gibbs_sampler(char **dna, int k, int t, int n, int motif_length, FILE* output);
double **profile_with_pseudocounts_without_i(char **motifs, int pseudocount, int motif_length, int number_of_motifs, int i_to_exclude, FILE* output);
int score(char **motifs, int motif_length, int number_of_motifs, FILE* output);
char *generate_random_kmer(char *motif_excluded, int k, int motif_length, int t, double **profile_matrix, FILE* output);
int *generate_unfair_dice_matrix(double *probabilities_matrix, int motif_length, int k);
int unfair_dice(double *probabilities_matrix, int *unfair_dice_matrix, int motif_length, int k);
double *generate_probabilities_of_all_kmers(char *motif_excluded, int k, int motif_length, int t, double **profile_matrix, FILE* output);
int lenHelper(unsigned x) ;

int main(int argc, char *argv[])
{
  FILE *input, *output;
  int k, t, n; // k is the size of k-mers ; n is the times we need to iterate ; t is the number of motifs;
  int i, j;
  int pseudocount = 0;
  char *first_motif_buffer, *check;
  char **motifs;
  int bytes_read, score_of_motifs;
  int motif_length, number_of_motifs;
  double **profile_matrix;

  // Get all the data from input file and store them in //
  // variables k = size of k-mers, t = number of motifs //
  // n = number of iterations and motifs = array that   //
  // contains all the motifs                            //

  input = fopen("input_brca1.txt", "r+");
  if (input == NULL)
  {
    printf("ERROR: Input not found!\n");
    exit(EXIT_FAILURE);
  }

  output = fopen("output.txt", "w+");
  if (output == NULL)
  {
    printf("ERROR: Output not found!\n");
    exit(EXIT_FAILURE);
  }

  fseek(input, 0 , SEEK_SET);
  bytes_read = 0;

  fscanf(input, "%d %d %d", &k, &t, &n);
  printf("k : %d, t : %d, N : %d\n",k, t, n);

  bytes_read = lenHelper(k) + 1 + lenHelper(t) + 1 + lenHelper(n) + 1;
  printf("Bytes read: %d\n", bytes_read);
  fseek(input, bytes_read , SEEK_SET);

  first_motif_buffer = (char*)malloc(5000*sizeof(char));
  check = fgets(first_motif_buffer, 5000*sizeof(char), input);
  if(check == NULL)
    exit(EXIT_FAILURE);

  motif_length = strlen(first_motif_buffer);
  first_motif_buffer[motif_length - 1] = '\0';
  motif_length = strlen(first_motif_buffer);

  printf("Motif length: %d\n", motif_length);
  //printf("first motif: %s\n", first_motif_buffer);

  number_of_motifs = 1;
  motifs = (char**)malloc(number_of_motifs * sizeof(char*));
  (*motifs) = (char*)malloc((motif_length + 1) * sizeof(char));
  strcpy((*motifs), first_motif_buffer);
  //printf("%s\n", (*motifs));
  free(first_motif_buffer);

  while (check != NULL)
  {
    number_of_motifs++;
    motifs = (char**)realloc(motifs, number_of_motifs * sizeof(char*));
    if(motifs==NULL)
    {
      printf("ERROR: In realloc\n");
    }
    (*(motifs+number_of_motifs - 1)) = (char*)malloc((motif_length + 1) * sizeof(char));

    check = fgets((*(motifs + number_of_motifs - 1)), (motif_length + 1)*sizeof(char), input);
    if(check == NULL)
    {
      number_of_motifs--;
      break;
    }

    *((*(motifs + number_of_motifs - 1)) + motif_length) = '\0';
    //printf("%s\n", (*(motifs + number_of_motifs - 1)));
    fseek(input, 1 , SEEK_CUR);
  }
  printf("Number of motifs: %d\n", number_of_motifs);

  /*profile_matrix = profile_with_pseudocounts_without_i(motifs, pseudocount, motif_length, number_of_motifs, 0);
  // print profile_matrix //
  for (i = 0 ; i < 4 ; i++)
  {
    for (j = 0; j < motif_length; j++)
    {
      printf("%.3lf ",*((*(profile_matrix + i)) + j));
    }
    printf("\n");
  }*/

  /*score_of_motifs = score(motifs, motif_length, number_of_motifs);
  printf("Score: %d\n", score_of_motifs);*/
 //free all

  gibbs_sampler(motifs, k, t, n, motif_length, output);

  fclose(input);
  fclose(output);
  return 0;
}

char **gibbs_sampler(char **dna, int k, int t, int n, int motif_length, FILE* output)
{
  int i, j, l, m, starting_position_of_kmer, i_to_exclude;
  char **motifs, **best_motifs;
  double **profile_matrix;
  char *most_probable_kmer;
  int score_of_motifs, score_of_best_motifs;
  clock_t start_t, end_t;
  char buffer[100];
  int written_in_buffer;

  srand(time(NULL));

  motifs = (char**)malloc(t*sizeof(char*));
  best_motifs = (char**)malloc(t*sizeof(char*));

  start_t = clock();
  for (i = 0 ; i < t ; i++)
  {
    *(motifs + i) = (char*)malloc((k+1)*sizeof(char*));

    starting_position_of_kmer = abs(rand() % (motif_length - k + 1));

    strncpy(*(motifs + i), (*(dna + i) + starting_position_of_kmer), k);
    *(*(motifs + i) + k) = '\0';
    //printf("Motifs: i: %d starting_position_of_kmer: %d %s\n", i, starting_position_of_kmer , *(motifs + i));

    *(best_motifs + i) = (char*)malloc((k+1)*sizeof(char*));

    strncpy(*(best_motifs + i), *(motifs + i), k);
    *(*(best_motifs + i) + k) = '\0';
    //printf("Best Motifs: i: %d starting_position_of_kmer: %d %s\n", i, starting_position_of_kmer , *(best_motifs + i));
  }
  end_t = clock();
  written_in_buffer = sprintf(buffer, "For#1: Function: gibbs_sampler, Total time: %lf\n", (double)(end_t - start_t)/ (double)CLOCKS_PER_SEC);
  //printf(" HEREEEE %s %d\n", buffer, written_in_buffer );
  fwrite(buffer, sizeof(char)* written_in_buffer, 1, output);

  printf("Motifs: \n");
  for (i = 0 ; i < t ; i++)
    printf("%s\n", *(motifs+i));
  /*printf("Best Motifs: \n");
  for (i = 0 ; i < t ; i++)
    printf("%s\n", *(best_motifs+i));*/

  start_t = clock();
  for (j = 0; j < n; j++)
  {
    i_to_exclude = abs(rand() % t);
    /*printf("Random motif: %d\n", i_to_exclude);
    printf("Motifs: \n");
    for (i = 0 ; i < t ; i++)
      printf("%s\n", *(motifs+i));*/
    profile_matrix = profile_with_pseudocounts_without_i(motifs, 1, k, t, i_to_exclude, output);

    // print profile matrix //
    /*for (m = 0 ; m < 4 ; m++)
    {
      for (l = 0; l < k; l++)
      {
        printf("%.3lf ",*((*(profile_matrix + m)) + l));
      }
      printf("\n");
    }
    printf("\n");*/

    most_probable_kmer = generate_random_kmer(*(dna + i_to_exclude), k, motif_length, t, profile_matrix, output);
    //printf("MOST PROBABLE KMER: %s\n",most_probable_kmer);
    strncpy(*(motifs + i_to_exclude), most_probable_kmer, k);
    *(*(motifs + i_to_exclude) + k) = '\0';

    score_of_motifs = score(motifs, k, t, output);
    score_of_best_motifs = score(best_motifs, k, t, output);
    //printf("M: %d, BM: %d\n", score_of_motifs, score_of_best_motifs);

    if(score_of_motifs < score_of_best_motifs)
    {
      //printf("YES!\n");
      strncpy(*(best_motifs + i_to_exclude), most_probable_kmer, k);
      *(*(best_motifs + i_to_exclude) + k) = '\0';

      /*printf("Motifs: \n");
      for (i = 0 ; i < t ; i++)
        printf("%s\n", *(best_motifs+i));*/
    }

    for (i = 0 ; i < 4 ; i++)
        free(*(profile_matrix + i));
    free(profile_matrix);
  }
  end_t = clock();
  written_in_buffer = sprintf(buffer, "For#2: Function: gibbs_sampler, Total time: %lf\n", (double)(end_t - start_t)/ (double)CLOCKS_PER_SEC);
  fwrite(buffer, sizeof(char)* written_in_buffer, 1, output);

  score_of_motifs = score(best_motifs, k, t, output);
  printf("Score of best motifss: %d\n", score_of_motifs);
  
  printf("Best Motifs: \n");
  for (i = 0 ; i < t ; i++)
  {
    printf("%s\n", *(best_motifs+i));
    free(*(best_motifs + i));
  }

  free(best_motifs);
  return best_motifs;
}

double **profile_with_pseudocounts_without_i(char **motifs, int pseudocount, int motif_length, int number_of_motifs, int i_to_exclude, FILE* output)
{
  int i, j, l;
  int count_a, count_c, count_g, count_t, sum;
  double **profile_matrix; // matrix of size 4 * motif_length ; each line corresponds to a single nucleotide    //
                          // in the order A, C, G, T ; each column corresponds to the number of appearances    //
                          // of a certain nucleotide (depends on the line) in the same column in motifs matrix //
                          // in this number we add the pseudocount to prevent zeros in the matrix and divide   //
                          // the result with 2 * number_of_motifs to obtain the corresponding probability      //
  clock_t start_t, end_t;
  char buffer[100];
  int written_in_buffer;

  // initialization and memory allocation of profile_matrix //

  profile_matrix = NULL;
  profile_matrix = (double**)malloc(4*sizeof(double*));

  start_t = clock();
  // parallelizable //
  for (i = 0 ; i < 4 ; i++)
  {
    *(profile_matrix + i) = NULL;
    *(profile_matrix + i) = (double*)malloc(motif_length*sizeof(double*));
    for (j = 0; j < motif_length; j++)
    {
      *((*(profile_matrix + i)) + j) = (double)pseudocount;
    }
  }
  end_t = clock();
  written_in_buffer = sprintf(buffer, "For#1: Function: profile_with_pseudocounts_without_i, Total time: %lf\n", (double)(end_t - start_t)/ (double)CLOCKS_PER_SEC);
  fwrite(buffer, sizeof(char)* written_in_buffer, 1, output);

  start_t = clock();
  // parallelizable //
  for (j = 0; j < motif_length; j++)
  {
    count_a = 0;
    count_c = 0;
    count_g = 0;
    count_t = 0;

    for(i = 0; i < number_of_motifs; i++)
    {
      if (i == i_to_exclude)
        continue;
      if((*((*(motifs+i))+j)) == 'A')
        count_a++;
      else if((*((*(motifs+i))+j)) == 'C')
        count_c++;
      else if((*((*(motifs+i))+j)) == 'G')
        count_g++;
      else if((*((*(motifs+i))+j)) == 'T')
        count_t++;
    }

    //printf("count_a = %d, count_c = %d, count_g = %d, count_t =%d\n", count_a, count_c, count_g, count_t);

    sum = count_a + count_c + count_g + count_t + 4*pseudocount;
    //printf("SUM: %d\n", sum);
    *((*(profile_matrix + 0)) + j) = *((*(profile_matrix + 0)) + j) + (double)count_a;
    //*((*(profile_matrix + 0)) + j) = *((*(profile_matrix + 0)) + j) / (double)sum;

    *((*(profile_matrix + 1)) + j) = *((*(profile_matrix + 1)) + j) + (double)count_c;
    //*((*(profile_matrix + 1)) + j) = *((*(profile_matrix + 1)) + j) / (double)sum;

    *((*(profile_matrix + 2)) + j) = *((*(profile_matrix + 2)) + j) + (double)count_g;
    //*((*(profile_matrix + 2)) + j) = *((*(profile_matrix + 2)) + j) / (double)sum;

    *((*(profile_matrix + 3)) + j) = *((*(profile_matrix + 3)) + j) + (double)count_t;
    //*((*(profile_matrix + 3)) + j) = *((*(profile_matrix + 3)) + j) / (double)sum;
  }
  end_t = clock();
  written_in_buffer = sprintf(buffer, "For#2: Function: profile_with_pseudocounts_without_i, Total time: %lf\n", (double)(end_t - start_t)/ (double)CLOCKS_PER_SEC);
  fwrite(buffer, sizeof(char)* written_in_buffer, 1, output);

  return profile_matrix;
}

int score(char **motifs, int motif_length, int number_of_motifs, FILE* output)
{
  int i, j;
  int count_a, count_c, count_g, count_t;
  int score; // the number of unpopular (lower case) letters in the motif matrix //
  clock_t start_t, end_t;
  char buffer[100];
  int written_in_buffer;

  score = 0;

  start_t = clock();
  // parallelizable //
  for (j = 0; j < motif_length; j++)
  {
    count_a = 0;
    count_c = 0;
    count_g = 0;
    count_t = 0;

    for(i = 0; i < number_of_motifs; i++)
    {
      if((*((*(motifs+i))+j)) == 'A')
        count_a++;
      else if((*((*(motifs+i))+j)) == 'C')
        count_c++;
      else if((*((*(motifs+i))+j)) == 'G')
        count_g++;
      else if((*((*(motifs+i))+j)) == 'T')
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

  end_t = clock();
  written_in_buffer = sprintf(buffer, "For#1: Function: score, Total time: %lf\n", (double)(end_t - start_t)/ (double)CLOCKS_PER_SEC);
  fwrite(buffer, sizeof(char)* written_in_buffer, 1, output);

  return score;
}

char *generate_random_kmer(char *motif_excluded, int k, int motif_length, int t, double **profile_matrix, FILE* output)
{
  double *probabilities_matrix;
  int most_probable_i, i;
  double most_probable_val;
  char *most_probable_kmer;
  int *unfair_dice_matrix;
  int unfair_dice_result;
  clock_t start_t, end_t;
  char buffer[100];
  int written_in_buffer;

  most_probable_kmer = (char*)malloc((k+1)*sizeof(char));

  probabilities_matrix = generate_probabilities_of_all_kmers(motif_excluded, k, motif_length, t, profile_matrix, output);

  most_probable_i = 0;
  most_probable_val = 0;

  start_t = clock();
  for (i = 0; i < motif_length - k + 1 ; i++)
  {
    if ((*(probabilities_matrix + i)) >= most_probable_val)
    {
      most_probable_val = (*(probabilities_matrix + i));
      most_probable_i = i;
    }
  }
  end_t = clock();
  written_in_buffer = sprintf(buffer, "For#1: Function: generate_random_kmer, Total time: %lf\n", (double)(end_t - start_t)/ (double)CLOCKS_PER_SEC);
  fwrite(buffer, sizeof(char)* written_in_buffer, 1, output);

  /*unfair_dice_matrix = generate_unfair_dice_matrix(probabilities_matrix, motif_length, k);
  unfair_dice_result = unfair_dice(probabilities_matrix, unfair_dice_matrix, motif_length, k);


  strncpy(most_probable_kmer, (motif_excluded + unfair_dice_result), k*sizeof(char));
  most_probable_kmer[k] = '\0';
*/
  strncpy(most_probable_kmer, (motif_excluded + most_probable_i), k*sizeof(char));
  most_probable_kmer[k] = '\0';
  //printf("MOST PROB KMER %s\n", most_probable_kmer);
  //free(unfair_dice_matrix);
  free(probabilities_matrix);
  //printf("MOST PROBABLE: %d %s %s, %lf\n",most_probable_i, motif_excluded, most_probable_kmer, most_probable_val);
  return most_probable_kmer;
}

int *generate_unfair_dice_matrix(double *probabilities_matrix, int motif_length, int k)
{
  unsigned long long int sum;
  int i, j;
  int *unfair_dice_matrix;

  sum = 0;
  for(i = 0; i < motif_length - k + 1; i++)
    sum = sum + (unsigned int)(*(probabilities_matrix + i));
  printf("%lld SUMM\n", sum);
  unfair_dice_matrix = (int*) malloc((abs(sum))*sizeof(int));
  for(i = 0; i < motif_length - k + 1; i++)
  {
    for(j = 0; j <  *(probabilities_matrix + i); j++)
      *(unfair_dice_matrix + j) = i;
  }
  return unfair_dice_matrix;
}

int unfair_dice(double *probabilities_matrix, int *unfair_dice_matrix, int motif_length, int k)
{
  int rand_pos, unfair_dice_result, sum, i;

  sum = 0;
  for(i = 0; i < motif_length - k + 1; i++)
    sum = sum + *(probabilities_matrix + i);

  rand_pos = abs(rand() % (sum - 1));
  printf("RAND POS = %d\n", rand_pos);
  unfair_dice_result = *(unfair_dice_matrix + rand_pos);
  printf("UNFAIR DICE RESULT = %d\n", unfair_dice_result);
  return unfair_dice_result;
}

double *generate_probabilities_of_all_kmers(char *motif_excluded, int k, int motif_length, int t, double **profile_matrix, FILE* output)
{
  int i, j;
  double *probabilities_matrix;
  clock_t start_t, end_t;
  char buffer[100];
  int written_in_buffer;

  probabilities_matrix = NULL;
  probabilities_matrix = (double*)malloc((motif_length - k + 1) * sizeof(double));

  start_t = clock();
  for (i = 0; i < motif_length - k + 1 ; i++)
  {
    *(probabilities_matrix + i) = 1;
    // parallelizable //
    for (j = 0; j < k; j++)
    {
      if (*(motif_excluded+i+j) == 'A')
        *(probabilities_matrix + i) = (*(probabilities_matrix + i)) * (*((*(profile_matrix + 0)) + j)) ;
      else if (*(motif_excluded+i+j) == 'C')
        *(probabilities_matrix + i) = (*(probabilities_matrix + i)) * (*((*(profile_matrix + 1)) + j));
      else if (*(motif_excluded+i+j) == 'G')
        *(probabilities_matrix + i) = (*(probabilities_matrix + i)) * (*((*(profile_matrix + 2)) + j));
      else if (*(motif_excluded+i+j) == 'T')
        *(probabilities_matrix + i) = (*(probabilities_matrix + i)) * (*((*(profile_matrix + 3)) + j));
      /*else
        printf("JUNK\n");*/
    }
    //*(probabilities_matrix + i) = *(probabilities_matrix + i) * (double)(1000000000);
  }
  // print probabilities_matrix for debuging //
  /*printf("Probabilities matrix: \n");
  for (i = 0; i < motif_length - k + 1 ; i++)
  {
    printf("%lf\n", *(probabilities_matrix + i));

  }*/
  end_t = clock();
  written_in_buffer = sprintf(buffer, "For#1: Function: generate_probabilities_of_all_kmers, Total time: %lf\n", (double)(end_t - start_t)/ (double)CLOCKS_PER_SEC);
  fwrite(buffer, sizeof(char)* written_in_buffer, 1, output);

  return probabilities_matrix;
}

int lenHelper(unsigned x)
{
    if (x >= 1000000000) return 10;
    if (x >= 100000000)  return 9;
    if (x >= 10000000)   return 8;
    if (x >= 1000000)    return 7;
    if (x >= 100000)     return 6;
    if (x >= 10000)      return 5;
    if (x >= 1000)       return 4;
    if (x >= 100)        return 3;
    if (x >= 10)         return 2;
    return 1;
}
