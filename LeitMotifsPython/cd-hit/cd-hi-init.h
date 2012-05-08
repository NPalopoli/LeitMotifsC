int NAA1 = MAX_UAA;
int NAA2 = NAA1 * NAA1;
int NAA3 = NAA1 * NAA1 * NAA1;
int NAA4 = NAA1 * NAA1 * NAA1 * NAA1;
int NAA5 = NAA1 * NAA1 * NAA1 * NAA1 * NAA1;
int NAA6 = NAA3 * NAA3;
int NAA7 = NAA3 * NAA4;
int NAA8 = NAA4 * NAA4;
int NAA9 = NAA4 * NAA5;
int NAA10= NAA5 * NAA5;
int NAA11= NAA5 * NAA6;
int NAA12= NAA6 * NAA6;
int NAAN_array[13] = {
  1, NAA1, NAA2, NAA3, NAA4, NAA5, NAA6, NAA7, NAA8, NAA9, NAA10, NAA11, NAA12
};
int NAA_top_limit = 5;

struct tms CPU_current, CPU_begin, CPU_end;


int DB_no, NR_no, NR90_no;
int *NR_len;
int *NR_idx;                // idx table -> unsorted no
int *NR90_idx;              // idx table -> old no
int *NR_clstr_no;           // cluster no.
int  *(*NR_info);           // other info, like alignment overlap, etc
char *NR_iden;
char *NR_seg;               // if SEG_no > 255, should be short
char *NR_flag;
char *(*NR_seq);
double NR_clstr = 0.9;
int NR_clstr100 = 90;
int NAA = 5;
int NAAN = NAA5;
int des_len = 20;
int tolerance = 2;
long long option_M  = 400000000;
int       option_l  = 10;
int       option_p  = 0;
int       option_g  = 0;
int       option_G  = 1;             
int       option_B  = 0;
int       option_b  = 20;
double    option_s  = 0.0;
int       option_S  = 99999999;
double    option_s2 = 1.0;
int       option_S2 = 0;
double    option_aL = 0.0;       
int       option_AL = 99999999;  
double    option_aS = 0.0;       
int       option_AS = 99999999;  
int       option_A  = 0;  
int       option_r  = 0;

IDX_TBL word_table;

int SEG_no;
int SEG_b[MAX_SEG], SEG_e[MAX_SEG], SEG90_b[MAX_SEG], SEG90_e[MAX_SEG];
char db_swap[MAX_SEG][MAX_FILE_NAME];

// following for mcd-hit
int *NR90f_idx;
int NR90f_no;
int NR_frag_no;
int Frag_size = 0;
int SEG90f_b[MAX_SEG], SEG90f_e[MAX_SEG];
// end

AA_MATRIX mat;

int naa_stat_start_percent = 40;
int naa_stat[5][61][4] = {

  // cover 0.99
  {
    // N=5   N=4   N=3   N=2
    {  0,    0,    0,    7,  },  // 40%
    {  0,    0,    0,    8,  },  // 41%
    {  0,    0,    0,    9,  },  // 42%
    {  0,    0,    0,    9,  },  // 43%
    {  0,    0,    1,   10,  },  // 44%
    {  0,    0,    1,   11,  },  // 45%
    {  0,    0,    1,   12,  },  // 46%
    {  0,    0,    2,   13,  },  // 47%
    {  0,    0,    2,   14,  },  // 48%
    {  0,    0,    4,   16,  },  // 49%
    {  0,    0,    4,   16,  },  // 50%
    {  0,    0,    5,   17,  },  // 51%
    {  0,    0,    5,   18,  },  // 52%
    {  0,    0,    7,   20,  },  // 53%
    {  0,    1,    7,   21,  },  // 54%
    {  0,    1,    7,   21,  },  // 55%
    {  0,    2,    8,   23,  },  // 56%
    {  0,    2,    8,   25,  },  // 57%
    {  0,    2,   10,   25,  },  // 58%
    {  0,    3,   10,   26,  },  // 59%
    {  0,    4,   13,   28,  },  // 60%
    {  0,    5,   13,   30,  },  // 61%
    {  0,    5,   14,   30,  },  // 62%
    {  1,    6,   15,   33,  },  // 63%
    {  2,    7,   17,   34,  },  // 64%
    {  2,    7,   17,   35,  },  // 65%
    {  2,    9,   20,   37,  },  // 66%
    {  4,   10,   20,   37,  },  // 67%
    {  4,   11,   22,   40,  },  // 68%
    {  5,   12,   24,   41,  },  // 69%
    {  5,   12,   25,   42,  },  // 70%
    {  6,   16,   27,   43,  },  // 71%
    {  8,   16,   27,   45,  },  // 72%
    {  9,   17,   29,   47,  },  // 73%
    { 10,   18,   31,   47,  },  // 74%
    { 10,   20,   32,   50,  },  // 75%
    { 12,   20,   32,   51,  },  // 76%
    { 14,   22,   36,   54,  },  // 77%
    { 15,   24,   37,   55,  },  // 78%
    { 17,   26,   41,   58,  },  // 79%
    { 18,   29,   41,   59,  },  // 80%
    { 20,   30,   45,   60,  },  // 81%
    { 24,   35,   48,   62,  },  // 82%
    { 26,   36,   48,   64,  },  // 83%
    { 27,   38,   51,   65,  },  // 84%
    { 31,   43,   54,   68,  },  // 85%
    { 35,   43,   55,   70,  },  // 86%
    { 36,   48,   60,   71,  },  // 87%
    { 36,   50,   61,   73,  },  // 88%
    { 40,   50,   61,   75,  },  // 89%
    { 45,   54,   65,   75,  },  // 90%
    { 52,   60,   70,   79,  },  // 91%
    { 53,   62,   71,   81,  },  // 92%
    { 57,   66,   75,   84,  },  // 93%
    { 57,   66,   76,   85,  },  // 94%
    { 64,   71,   78,   85,  },  // 95%
    { 70,   75,   82,   89,  },  // 96%
    { 77,   81,   86,   92,  },  // 97%
    { 82,   86,   90,   94,  },  // 98%
    { 83,   87,   91,   95,  },  // 99%
    { 91,   93,   95,   97,  },  // 100%
  },
  // cover 0.95
  {
    // N=5   N=4   N=3   N=2
    {  0,    0,    1,    9,  },  // 40%
    {  0,    0,    2,   10,  },  // 41%
    {  0,    0,    2,   11,  },  // 42%
    {  0,    0,    3,   12,  },  // 43%
    {  0,    0,    3,   12,  },  // 44%
    {  0,    0,    4,   14,  },  // 45%
    {  0,    0,    4,   14,  },  // 46%
    {  0,    1,    5,   16,  },  // 47%
    {  0,    1,    6,   17,  },  // 48%
    {  0,    2,    7,   19,  },  // 49%
    {  0,    2,    8,   19,  },  // 50%
    {  0,    2,    8,   20,  },  // 51%
    {  0,    2,    9,   21,  },  // 52%
    {  0,    4,   10,   23,  },  // 53%
    {  1,    4,   11,   24,  },  // 54%
    {  1,    4,   11,   24,  },  // 55%
    {  1,    5,   13,   26,  },  // 56%
    {  2,    5,   13,   27,  },  // 57%
    {  2,    6,   15,   29,  },  // 58%
    {  2,    7,   15,   30,  },  // 59%
    {  3,    8,   16,   31,  },  // 60%
    {  4,    8,   18,   32,  },  // 61%
    {  4,    9,   18,   33,  },  // 62%
    {  5,   11,   20,   36,  },  // 63%
    {  6,   12,   22,   37,  },  // 64%
    {  6,   12,   22,   38,  },  // 65%
    {  8,   14,   24,   40,  },  // 66%
    {  8,   15,   25,   41,  },  // 67%
    { 10,   16,   27,   42,  },  // 68%
    { 10,   18,   28,   45,  },  // 69%
    { 11,   18,   29,   45,  },  // 70%
    { 14,   21,   31,   47,  },  // 71%
    { 14,   22,   32,   48,  },  // 72%
    { 14,   22,   33,   50,  },  // 73%
    { 17,   24,   36,   52,  },  // 74%
    { 17,   25,   36,   52,  },  // 75%
    { 18,   27,   39,   54,  },  // 76%
    { 20,   29,   41,   56,  },  // 77%
    { 21,   31,   42,   58,  },  // 78%
    { 21,   31,   46,   60,  },  // 79%
    { 27,   35,   46,   60,  },  // 80%
    { 28,   37,   50,   63,  },  // 81%
    { 31,   38,   50,   64,  },  // 82%
    { 34,   43,   53,   66,  },  // 83%
    { 36,   45,   54,   67,  },  // 84%
    { 41,   50,   60,   70,  },  // 85%
    { 43,   51,   60,   71,  },  // 86%
    { 45,   54,   63,   74,  },  // 87%
    { 48,   55,   64,   75,  },  // 88%
    { 54,   60,   68,   78,  },  // 89%
    { 55,   62,   71,   80,  },  // 90%
    { 56,   63,   71,   80,  },  // 91%
    { 64,   70,   76,   84,  },  // 92%
    { 69,   74,   80,   86,  },  // 93%
    { 73,   78,   83,   88,  },  // 94%
    { 74,   78,   84,   89,  },  // 95%
    { 80,   84,   87,   91,  },  // 96%
    { 83,   86,   90,   93,  },  // 97%
    { 86,   89,   92,   95,  },  // 98%
    { 91,   93,   95,   97,  },  // 99%
    { 92,   93,   95,   97,  },  // 100%
  },
  // cover 0.9
  {
    // N=5   N=4   N=3   N=2
    {  0,    0,    2,   11,  },  // 40%
    {  0,    0,    3,   12,  },  // 41%
    {  0,    0,    3,   12,  },  // 42%
    {  0,    1,    4,   13,  },  // 43%
    {  0,    1,    5,   14,  },  // 44%
    {  0,    1,    5,   15,  },  // 45%
    {  0,    1,    6,   16,  },  // 46%
    {  0,    2,    7,   18,  },  // 47%
    {  0,    2,    7,   18,  },  // 48%
    {  0,    3,    9,   20,  },  // 49%
    {  1,    4,    9,   20,  },  // 50%
    {  1,    4,   10,   21,  },  // 51%
    {  1,    4,   11,   23,  },  // 52%
    {  2,    5,   12,   24,  },  // 53%
    {  2,    5,   12,   25,  },  // 54%
    {  2,    6,   13,   26,  },  // 55%
    {  3,    7,   14,   28,  },  // 56%
    {  3,    7,   15,   28,  },  // 57%
    {  4,    8,   16,   30,  },  // 58%
    {  5,    9,   17,   31,  },  // 59%
    {  5,   10,   18,   32,  },  // 60%
    {  6,   11,   20,   35,  },  // 61%
    {  6,   11,   20,   35,  },  // 62%
    {  7,   13,   22,   38,  },  // 63%
    {  8,   14,   23,   39,  },  // 64%
    {  8,   15,   24,   39,  },  // 65%
    { 10,   16,   26,   42,  },  // 66%
    { 10,   17,   27,   42,  },  // 67%
    { 12,   19,   29,   44,  },  // 68%
    { 13,   20,   30,   46,  },  // 69%
    { 13,   21,   31,   47,  },  // 70%
    { 16,   23,   33,   48,  },  // 71%
    { 18,   25,   34,   50,  },  // 72%
    { 18,   26,   36,   51,  },  // 73%
    { 19,   28,   38,   53,  },  // 74%
    { 20,   29,   38,   53,  },  // 75%
    { 23,   30,   41,   56,  },  // 76%
    { 24,   33,   43,   57,  },  // 77%
    { 26,   34,   45,   59,  },  // 78%
    { 28,   37,   48,   61,  },  // 79%
    { 30,   37,   48,   62,  },  // 80%
    { 33,   42,   52,   64,  },  // 81%
    { 35,   43,   53,   65,  },  // 82%
    { 38,   47,   56,   68,  },  // 83%
    { 40,   47,   56,   68,  },  // 84%
    { 44,   53,   61,   71,  },  // 85%
    { 45,   53,   62,   73,  },  // 86%
    { 50,   58,   66,   75,  },  // 87%
    { 51,   58,   66,   76,  },  // 88%
    { 57,   63,   71,   79,  },  // 89%
    { 60,   66,   72,   81,  },  // 90%
    { 62,   68,   75,   83,  },  // 91%
    { 70,   74,   80,   85,  },  // 92%
    { 74,   78,   82,   88,  },  // 93%
    { 85,   87,   90,   92,  },  // 94%
    { 86,   88,   90,   92,  },  // 95%
    { 87,   89,   91,   93,  },  // 96%
    { 87,   89,   92,   94,  },  // 97%
    { 89,   91,   93,   96,  },  // 98%
    { 93,   94,   96,   97,  },  // 99%
    { 94,   95,   97,   98,  },  // 100%
  },
  // cover 0.8
  {
    // N=5   N=4   N=3   N=2
    {  0,    1,    4,   13,  },  // 40%
    {  0,    1,    5,   13,  },  // 41%
    {  0,    1,    5,   14,  },  // 42%
    {  0,    2,    6,   15,  },  // 43%
    {  0,    2,    6,   16,  },  // 44%
    {  0,    2,    7,   17,  },  // 45%
    {  1,    3,    8,   18,  },  // 46%
    {  1,    4,    9,   20,  },  // 47%
    {  1,    4,    9,   20,  },  // 48%
    {  2,    5,   11,   22,  },  // 49%
    {  2,    5,   11,   22,  },  // 50%
    {  2,    6,   12,   24,  },  // 51%
    {  3,    6,   13,   25,  },  // 52%
    {  3,    7,   14,   26,  },  // 53%
    {  4,    8,   14,   27,  },  // 54%
    {  4,    8,   15,   28,  },  // 55%
    {  5,    9,   17,   30,  },  // 56%
    {  5,    9,   17,   30,  },  // 57%
    {  6,   11,   19,   32,  },  // 58%
    {  7,   12,   20,   34,  },  // 59%
    {  8,   12,   20,   34,  },  // 60%
    {  9,   14,   22,   37,  },  // 61%
    {  9,   14,   23,   37,  },  // 62%
    { 10,   16,   25,   39,  },  // 63%
    { 11,   17,   26,   41,  },  // 64%
    { 12,   18,   27,   41,  },  // 65%
    { 13,   20,   28,   43,  },  // 66%
    { 14,   21,   30,   45,  },  // 67%
    { 15,   22,   31,   46,  },  // 68%
    { 17,   24,   33,   48,  },  // 69%
    { 17,   24,   34,   48,  },  // 70%
    { 19,   26,   36,   50,  },  // 71%
    { 20,   27,   37,   51,  },  // 72%
    { 21,   29,   39,   53,  },  // 73%
    { 23,   31,   41,   55,  },  // 74%
    { 23,   31,   41,   55,  },  // 75%
    { 26,   34,   44,   58,  },  // 76%
    { 28,   36,   46,   59,  },  // 77%
    { 29,   37,   47,   60,  },  // 78%
    { 34,   41,   50,   62,  },  // 79%
    { 34,   42,   51,   63,  },  // 80%
    { 38,   45,   55,   66,  },  // 81%
    { 39,   46,   55,   67,  },  // 82%
    { 44,   51,   60,   70,  },  // 83%
    { 44,   51,   60,   70,  },  // 84%
    { 49,   56,   64,   73,  },  // 85%
    { 50,   57,   64,   74,  },  // 86%
    { 57,   63,   69,   77,  },  // 87%
    { 58,   64,   70,   78,  },  // 88%
    { 68,   71,   76,   82,  },  // 89%
    { 68,   72,   77,   83,  },  // 90%
    { 75,   79,   81,   85,  },  // 91%
    { 86,   87,   89,   90,  },  // 92%
    { 88,   89,   90,   92,  },  // 93%
    { 90,   91,   92,   93,  },  // 94%
    { 91,   92,   93,   94,  },  // 95%
    { 92,   94,   94,   95,  },  // 96%
    { 93,   94,   95,   96,  },  // 97%
    { 94,   95,   95,   96,  },  // 98%
    { 94,   95,   96,   98,  },  // 99%
    { 95,   96,   97,   98,  },  // 100%
  },
  // cover 0.6
  {
    // N=5   N=4   N=3   N=2
    {  1,    2,    6,   15,  },  // 40%
    {  1,    3,    7,   16,  },  // 41%
    {  1,    3,    8,   17,  },  // 42%
    {  2,    4,    9,   18,  },  // 43%
    {  2,    4,    9,   19,  },  // 44%
    {  2,    5,   10,   20,  },  // 45%
    {  3,    5,   10,   21,  },  // 46%
    {  3,    6,   12,   22,  },  // 47%
    {  3,    6,   12,   23,  },  // 48%
    {  4,    8,   14,   25,  },  // 49%
    {  4,    8,   14,   25,  },  // 50%
    {  5,    8,   15,   26,  },  // 51%
    {  5,    9,   16,   27,  },  // 52%
    {  6,   10,   17,   29,  },  // 53%
    {  6,   11,   18,   30,  },  // 54%
    {  7,   11,   18,   31,  },  // 55%
    {  8,   12,   20,   32,  },  // 56%
    {  8,   13,   20,   33,  },  // 57%
    { 10,   14,   22,   35,  },  // 58%
    { 10,   15,   23,   37,  },  // 59%
    { 11,   16,   24,   37,  },  // 60%
    { 12,   18,   26,   39,  },  // 61%
    { 13,   18,   26,   40,  },  // 62%
    { 14,   20,   28,   42,  },  // 63%
    { 16,   22,   30,   43,  },  // 64%
    { 16,   22,   31,   44,  },  // 65%
    { 17,   23,   32,   45,  },  // 66%
    { 18,   25,   33,   47,  },  // 67%
    { 19,   26,   35,   48,  },  // 68%
    { 21,   27,   36,   50,  },  // 69%
    { 22,   29,   37,   51,  },  // 70%
    { 24,   30,   39,   52,  },  // 71%
    { 25,   32,   41,   53,  },  // 72%
    { 26,   33,   42,   55,  },  // 73%
    { 29,   35,   44,   57,  },  // 74%
    { 29,   36,   45,   57,  },  // 75%
    { 32,   39,   48,   60,  },  // 76%
    { 34,   41,   50,   61,  },  // 77%
    { 36,   43,   51,   62,  },  // 78%
    { 40,   46,   54,   65,  },  // 79%
    { 40,   46,   54,   65,  },  // 80%
    { 46,   52,   59,   68,  },  // 81%
    { 46,   52,   60,   69,  },  // 82%
    { 53,   59,   65,   73,  },  // 83%
    { 54,   60,   66,   73,  },  // 84%
    { 63,   67,   73,   78,  },  // 85%
    { 68,   71,   75,   79,  },  // 86%
    { 78,   80,   82,   85,  },  // 87%
    { 79,   81,   83,   85,  },  // 88%
    { 83,   85,   86,   87,  },  // 89%
    { 85,   86,   87,   89,  },  // 90%
    { 86,   88,   89,   90,  },  // 91%
    { 88,   89,   90,   91,  },  // 92%
    { 90,   90,   91,   92,  },  // 93%
    { 91,   92,   92,   93,  },  // 94%
    { 92,   93,   94,   94,  },  // 95%
    { 94,   94,   95,   95,  },  // 96%
    { 95,   95,   96,   96,  },  // 97%
    { 95,   96,   97,   97,  },  // 98%
    { 96,   96,   97,   98,  },  // 99%
    { 97,   98,   98,   99,  },  // 100%
  },
};
