// =============================================================================
// CD-HI-EST
// http://cd-hit.org/
// Cluster Database at High Identity (EST version)
// modified from CD-HI
//
// program written by 
//                                      Weizhong Li
//                                      UCSD, San Diego Supercomputer Center
//                                      La Jolla, CA, 92093
//                                      Email liwz@sdsc.edu
//                 at
//                                      Adam Godzik's lab
//                                      The Burnham Institute
//                                      La Jolla, CA, 92037
//                                      Email adam@burnham-inst.org
// =============================================================================

#include "cd-hi.h"
//over-write some defs in cd-hi.h for est version
#undef MAX_UAA
#define MAX_UAA 4

#include "cd-hi-init.h"
//over-write some defs in cd-hi-init.h for est version
int *Comp_AAN_idx;
char seqi_comp[MAX_SEQ];
int aan_list_comp[MAX_SEQ];

void setaa_to_na();
void make_comp_short_word_index(int NAA, int *NAAN_array, int *Comp_AAN_idx);
void make_comp_iseq(int len, char *iseq_comp, char *iseq);


int DB2_no;
int NR2_no;
int *NR2_len;
int *NR2_clstr_no;
int NR2_red_no;
int *(*NR2_info);
char *NR2_iden;
char *NR2_flag;
char *(*NR2_seq);
char *NR2_seg;
int SEG2_no;
int SEG2_b[MAX_SEG], SEG2_e[MAX_SEG];
// next two control how if seqs in db2 is longer than reps in db1
// by deault, only seqs in db2 that are shorter than rep in db1
// are clustered to the rep in db1


////////////////////////////////////  MAIN /////////////////////////////////////
int main(int argc, char **argv) {
  int i, j, k, i1, j1, k1, i0, j0, k0, sg_i, sg_j;
  int si, sj, sk;
  char db_in[MAX_FILE_NAME];
  char db2_in[MAX_FILE_NAME];
  char db_out[MAX_FILE_NAME];
  char db_clstr[MAX_FILE_NAME];
  char db_clstr_bak[MAX_FILE_NAME];
  char db_clstr_old[MAX_FILE_NAME];
  char db_bin_swap[MAX_FILE_NAME];
  char db2_bin_swap[MAX_FILE_NAME];

  NAA = 8;
  NAAN = NAA8;
  NAA_top_limit = 12;
  setaa_to_na();
  mat.set_to_na(); // mat.set_gap(-6,-1);

  times(&CPU_begin);

  // ***********************************    parse command line and open file
  if (argc < 7) print_usage_est_2d(argv[0]);
  for (i=1; i<argc; i++) {
    if      (strcmp(argv[i], "-i" ) == 0) strncpy(db_in,  argv[++i], MAX_FILE_NAME-1);
    else if (strcmp(argv[i], "-i2") == 0) strncpy(db2_in, argv[++i], MAX_FILE_NAME-1);
    else if (strcmp(argv[i], "-o" ) == 0) strncpy(db_out, argv[++i], MAX_FILE_NAME-1);
    else if (strcmp(argv[i], "-M" ) == 0) option_M  = atoll(argv[++i]) * 1000000;
    else if (strcmp(argv[i], "-l" ) == 0) option_l  = atoi(argv[++i]);
    else if (strcmp(argv[i], "-c" ) == 0) NR_clstr  = atof(argv[++i]);
    else if (strcmp(argv[i], "-b" ) == 0) option_b  = atoi(argv[++i]);
    else if (strcmp(argv[i], "-n" ) == 0) NAA       = atoi(argv[++i]);
    else if (strcmp(argv[i], "-d" ) == 0) des_len   = atoi(argv[++i]);
    else if (strcmp(argv[i], "-s" ) == 0) option_s  = atof(argv[++i]);
    else if (strcmp(argv[i], "-S" ) == 0) option_S  = atoi(argv[++i]);
    else if (strcmp(argv[i], "-s2") == 0) option_s2 = atof(argv[++i]); 
    else if (strcmp(argv[i], "-S2") == 0) option_S2 = atoi(argv[++i]); 
    else if (strcmp(argv[i], "-B" ) == 0) option_B  = atoi(argv[++i]);
    else if (strcmp(argv[i], "-r")  == 0) option_r  = atoi(argv[++i]);
    else if (strcmp(argv[i], "-p" ) == 0) option_p  = atoi(argv[++i]);
    else if (strcmp(argv[i], "-g" ) == 0) option_g  = atoi(argv[++i]);
    else if (strcmp(argv[i], "-G" ) == 0) option_G  = atoi(argv[++i]);
    else if (strcmp(argv[i], "-aL") == 0) option_aL = atof(argv[++i]);
    else if (strcmp(argv[i], "-AL") == 0) option_AL = atoi(argv[++i]);
    else if (strcmp(argv[i], "-aS") == 0) option_aS = atof(argv[++i]);
    else if (strcmp(argv[i], "-AS") == 0) option_AS = atoi(argv[++i]);
    else if (strcmp(argv[i], "-A" ) == 0) option_A  = atoi(argv[++i]);
    else                                  print_usage_est_2d(argv[0]);
  }
  if (1) {
      if ((NR_clstr > 1.0) || (NR_clstr < 0.8)) 
        bomb_error("invalid clstr threshold, should >=0.8");
      if (option_b < 1 ) bomb_error("invalid band width");
      if (NAA < 2 || NAA > NAA_top_limit) bomb_error("invalid word length");
      if (des_len < 0) 
        bomb_error("too short description, not enough to identify sequences");
      if ((option_s<0) || (option_s>1)) 
        bomb_error("invalid value for -s");
      if (option_S<0) bomb_error("invalid value for -S");
      if ((option_s2<0) || (option_s2>1))
        bomb_error("invalid value for -s");
      if (option_S2<0) bomb_error("invalid value for -S");
      if (option_G == 0) option_p = 1;
      if (option_aS < option_aL) option_aS = option_aL;
      if (option_AS > option_AL) option_AS = option_AL;
      if ((option_G == 0) && (option_aS == 0.0) && (option_A == 0))
        bomb_error("You are using local identity, but no -aS -aL -A option");
  }
  NR_clstr100 = (int) (NR_clstr * 100 );

  db_clstr[0]=0; strcat(db_clstr,db_out); strcat(db_clstr,".clstr");
  db_clstr_bak[0]=0;
  strcat(db_clstr_bak,db_out); strcat(db_clstr_bak,".bak.clstr");

  NAAN = NAAN_array[NAA];
  word_table.set_dna();
  word_table.init(NAA, NAAN);

  if (1) {
    if      ( NR_clstr > 0.9  && NAA < 8)
      cout << "Your word length is " << NAA
           << ", using 8 may be faster!" <<endl;
    else if ( NR_clstr > 0.87 && NAA < 5)
      cout << "Your word length is " << NAA
           << ", using 5 may be faster!" <<endl;
    else if ( NR_clstr > 0.80 && NAA < 4 )
      cout << "Your word length is " << NAA
           << ", using 4 may be faster!" <<endl;
    else if ( NR_clstr > 0.75 && NAA < 3 )
      cout << "Your word length is " << NAA
           << ", using 3 may be faster!" <<endl;
  }

  if ( option_l <= NAA ) bomb_error("Too short -l, redefine it");

  if ( option_r ) {
    if ((Comp_AAN_idx = new int[NAAN]) == NULL) bomb_error("Memory");
    make_comp_short_word_index(NAA, NAAN_array, Comp_AAN_idx);
  }


  ifstream in1a(db_in);    if (! in1a) bomb_error("Can not open", db_in);
  ifstream in2a(db2_in);   if (! in2a) bomb_error("Can not open ", db2_in);
  ofstream out1(db_out);   if (! out1) bomb_error("Can not open", db_out);
  ofstream out2(db_clstr); if (! out2) bomb_error("Can not open", db_clstr);
  ofstream out2b(db_clstr_bak); if (! out2b)
                                  bomb_error("Can not open", db_clstr_bak);
  strcpy(db_bin_swap,db_out); strcat(db_bin_swap,".BINSWAP");
  strcpy(db2_bin_swap,db_out); strcat(db2_bin_swap,".2.BINSWAP");

  DB_no = db_seq_no_test(in1a); in1a.close();
  DB2_no= db_seq_no_test(in2a); in2a.close();

  ifstream in1(db_in);
  if (! in1) bomb_error("Can not open file", db_in); 
  ifstream in2(db2_in);
  if (! in2) bomb_error("Can not open file", db2_in);

  if ((NR_len      = new int   [DB_no]) == NULL) bomb_error("Memory");
  if ((NR_idx      = new int   [DB_no]) == NULL) bomb_error("Memory");
  if ((NR90_idx    = new int   [DB_no]) == NULL) bomb_error("Memory");
  if ((NR_clstr_no = new int   [DB_no]) == NULL) bomb_error("Memory");
  if ((NR_iden     = new char  [DB_no]) == NULL) bomb_error("Memory");
  if ((NR_seg      = new char  [DB_no]) == NULL) bomb_error("Memory");
  if ((NR_flag     = new char  [DB_no]) == NULL) bomb_error("Memory");
  if ((NR_seq      = new char *[DB_no]) == NULL) bomb_error("Memory");
  int *Clstr_no, *Clstr_no_db1, *(*Clstr_list);
  if ((Clstr_no    = new int   [DB_no]) == NULL) bomb_error("Memory");
  if ((Clstr_no_db1= new int   [DB_no]) == NULL) bomb_error("Memory");
  if ((Clstr_list  = new int  *[DB_no]) == NULL) bomb_error("Memory");

  if ((NR2_len     = new int  [DB2_no]) == NULL) bomb_error("Memory");
  if ((NR2_clstr_no= new int  [DB2_no]) == NULL) bomb_error("Memory");
  if ((NR2_iden    = new char [DB2_no]) == NULL) bomb_error("Memory");
  if ((NR2_seg     = new char [DB2_no]) == NULL) bomb_error("Memory");
  if ((NR2_flag    = new char [DB2_no]) == NULL) bomb_error("Memory");
  if ((NR2_seq     = new char*[DB2_no]) == NULL) bomb_error("Memory");
  if (option_p) {
    if ((NR2_info  = new int *[DB2_no]) == NULL) bomb_error("Memory");
  }

  db_read_in(in1, db_bin_swap, option_B, option_l, NR_no, NR_seq,
             NR_len);
  in1.close();
  cout << "total seq in db1: " << NR_no << endl;
                                                                                   
  db_read_in(in2, db2_bin_swap, option_B, option_l, NR2_no, NR2_seq,
             NR2_len);
  in2.close();
  cout << "total seq in db2: " << NR2_no << endl;


  // ********************************************* init NR_flag
  for (i=0; i<NR_no; i++) NR_flag[i] = 0;
  for (i=0; i<DB2_no;i++) { NR2_flag[i]=0; NR2_iden[i]=0;
    NR2_clstr_no[i]=999999999;
  }

  sort_seqs_divide_segs(option_B, NR_no, NR_len, NR_idx, NR_seq, option_M,
                        NAAN, SEG_no, SEG_b, SEG_e, db_swap,db_out);

  for (sg_i=0; sg_i<SEG_no; sg_i++) {
    for (i1=SEG_b[sg_i]; i1<=SEG_e[sg_i]; i1++) {
      i = NR_idx[i1];
      NR_seg[i] = sg_i;
    }
  }

  db2_seqs_divide_segs(option_B, NR2_no, NR2_len, NR2_seq, option_M,
                       NAAN, SEG2_no, SEG2_b, SEG2_e);
  for (sg_i=0; sg_i<SEG2_no; sg_i++) {
    for (i1=SEG2_b[sg_i]; i1<=SEG2_e[sg_i]; i1++) {
      i = i1;
      NR2_seg[i] = sg_i;
    }
  }

  // *********************************************                Main loop
  char *seqi;
  double aa1_cutoff = NR_clstr;
  double aas_cutoff = 1 - (1-NR_clstr)*4;
  double aan_cutoff = 1 - (1-NR_clstr)*NAA;
  int len, hit_no, has_aas, iden_no, aan_no, segb, lens, len_tmp;
  int aan_list[MAX_SEQ];
  INTs aan_list_no[MAX_SEQ];
  INTs *look_and_count;
  int len_upper_bound, len_lower_bound;
  int check_this_info[32];
  double check_this_infod[32];
  check_this_info[0]   = option_p;
  check_this_info[10]  = option_G;
  check_this_infod[11] = option_aL;
  check_this_info[12]  = option_AL;
  check_this_infod[13] = option_aS;
  check_this_info[14]  = option_AS;
  check_this_info[15]  = option_A;
  check_this_info[20]  = option_g;
  if ((look_and_count= new INTs[NR_no]) == NULL) bomb_error("Memory");

  // write index table for database1
  cout << "compute index table for first database" << endl;
  NR90_no = 0;
  for (sg_i=0; sg_i<SEG_no; sg_i++) {
    if (SEG_no >1)
      cout << "SEG " << sg_i << " " << SEG_b[sg_i] << " " << SEG_e[sg_i] <<endl;
    if(option_B) read_swap_iseq(sg_i, db_bin_swap, NR_no, NR_seg, NR_seq);
                                                                                   
    word_table.clean();
    segb = NR90_no;
    for (i1=SEG_b[sg_i]; i1<=SEG_e[sg_i]; i1++) {
      i = NR_idx[i1];
      len = NR_len[i]; seqi = NR_seq[i];
      calc_ann_list(len, seqi, NAA, aan_no, aan_list, aan_list_no);
                                                                                   
      NR90_idx[NR90_no] = i;
      NR_clstr_no[i] = NR90_no; // positive value for representatives
      NR_iden[i] = 0;
      NR_flag[i] |= IS_REP;
      word_table.add_word_list2(aan_no, aan_list, aan_list_no, NR90_no);
      NR90_no++;
                                                                                   
      if ( (i1+1) % 100 == 0 ) {
        cout << ".";
        if ( (i1+1) % 1000 == 0 ) cout << i1+1 << " finished\t" << endl;
      }
    } // for (i1=SEG_b[sg_i]; i1<=SEG_e[sg_i]; i1++) {
                                                                                   
    SEG90_b[sg_i] = segb;  SEG90_e[sg_i] = NR90_no-1;
    if (SEG_no>1) word_table.write_tbl( db_swap[sg_i] );
    if(option_B) free_swap_iseq(sg_i, NR_no, NR_seg, NR_seq);
  } // for (sg_i=0; sg_i<SEG_no; sg_i++) {


  // compare database2 to database1
  NR2_red_no = 0;
  for (sg_i=0; sg_i<SEG2_no; sg_i++) {
    if (SEG2_no >1)
      cout << "SEG " << sg_i << " " << SEG2_b[sg_i] << " "
           << SEG2_e[sg_i] <<endl;
    if(option_B) read_swap_iseq(sg_i, db2_bin_swap, NR2_no, NR2_seg, NR2_seq);

    for (sg_j=0; sg_j<SEG_no; sg_j++) {
      cout << "Reading swap" << endl;
      // reading old segment
      if(option_B) read_swap_iseq(sg_j, db_bin_swap, NR_no, NR_seg, NR_seq);
      if (SEG_no>1) word_table.read_tbl(db_swap[sg_j]);
      cout << "Comparing with SEG " << sg_j << endl;

      for (i1=SEG2_b[sg_i]; i1<=SEG2_e[sg_i]; i1++) {
        i = i1;
        if (NR2_flag[i] & IS_REDUNDANT  ) continue;

        len = NR2_len[i]; seqi = NR2_seq[i];
        len_upper_bound = upper_bound_length_rep(len,
          option_s, option_S, option_aL, option_AL);

        len_lower_bound = len - option_S2;
        len_tmp = (int) ( ((double)len) * option_s2);
        if (len_tmp < len_lower_bound) len_lower_bound = len_tmp;
        has_aas = 0;
        iden_no =  NR2_iden[i];
        int flag = check_this_2d(len, seqi, has_aas,
               NAA, aan_no, aan_list, aan_list_no, look_and_count,
               hit_no, SEG90_b[sg_j], SEG90_e[sg_j], iden_no,
               aa1_cutoff, aas_cutoff, aan_cutoff,
               NR2_flag[i], NR_flag, len_upper_bound,
               lens, len_lower_bound, check_this_info, check_this_infod);

        if ((flag == 1) || (flag == -1)) {
          if (! option_g) {
            if (! option_B) delete [] NR2_seq[i];
            NR2_flag[i] |= IS_REDUNDANT ;
            NR2_red_no++;
          }
          NR2_clstr_no[i] = -hit_no-1;  // (-hit_no-1) for non representatives
          NR2_iden[i] = iden_no;
          if (flag == -1) NR2_flag[i] |= IS_MINUS_STRAND; // minus strand
          if (option_p){
            if ((NR2_info[i] = new int [4]) == NULL) bomb_error("Memory");
            NR2_info[i][0] = check_this_info[1]+1;
            NR2_info[i][1] = check_this_info[2]+1;
            NR2_info[i][2] = check_this_info[3]+1;
            NR2_info[i][3] = check_this_info[4]+1;
          }
        }

        if ( (i1+1) % 100 == 0 ) {
          cout << ".";
          if ( (i1+1) % 1000 == 0 )
            cout << i1+1 << " compared\t" << NR2_red_no << " clustered" << endl;
        }
      } //for (i1=SEG2_b[sg_i]; i1<=SEG2_e[sg_i]; i1++)
      if(option_B) free_swap_iseq(sg_j, NR_no, NR_seg, NR_seq);
    } // for (sg_j=SEG_no-1; sg_j>=0; sg_j--) {
                                                                                   
    if(option_B) free_swap_iseq(sg_i, NR2_no, NR2_seg, NR2_seq);
  } // for (sg_i=0; sg_i<SEG2_no; sg_i++) {

  if (option_g) {//delete redundant sequences in option_g mode
    for (i=0; i<NR2_no; i++)
      if (NR2_iden[i] > 0) {
        if (! option_B) delete [] NR2_seq[i];
        NR2_flag[i] |= IS_REDUNDANT ;
        NR2_red_no++;
      }
  }

  cout << endl;
  cout << NR2_no << " compared\t" << NR2_red_no << " clustered" << endl;

  if (! option_B) for (i=0; i<NR90_no; i++)  delete [] NR_seq[ NR90_idx[i] ]; 
  if (! option_B) {
    for (i=0; i<NR2_no; i++) {
      if (! (NR2_flag[i] & IS_REDUNDANT )) delete [] NR2_seq[i];
    }
  }

  cout << "writing non-redundant sequences from db2" << endl;
  ifstream in2b(db2_in);
  if ( ! in2b) bomb_error("Can not open file twice",db2_in);
  db_read_and_write(in2b, out1,option_l, des_len, NR2_seq, NR2_clstr_no);
  in2b.close(); out1.close();

  ifstream in1b(db_in);
  if ( ! in1b) bomb_error("Can not open file twice",db_in);
  db_read_des(in1b, option_l, des_len, NR_seq);
  in1b.close();


  // write a backup clstr file in case next step crashes
  for (i=0; i<NR_no; i++) {
    j1 = NR_clstr_no[i];
    if ( j1 < 0 ) j1 =-j1-1;
    out2b << j1 << "\t" << NR_len[i] << "aa, "<< NR_seq[i] << "...";
    if ( NR_iden[i]>0 ) out2b << " at " << int(NR_iden[i]) << "%" << endl;
    else                out2b << " *" << endl;
  }
  for (i=0; i<NR2_no; i++) {
    j1 = NR2_clstr_no[i];
    if ( j1 >= 0 ) continue; // skip nr seq from db2
    if ( j1 < 0 ) j1 =-j1-1;
    out2b << j1 << "\t" << NR2_len[i] << "aa, "<< NR2_seq[i] << "...";
    if ( NR2_iden[i]>0 ) {
      out2b << " at ";
      if (option_p)
        out2b << NR2_info[i][0] << ":" << NR2_info[i][1] << ":"
              << NR2_info[i][2] << ":" << NR2_info[i][3] << "/";
      out2b << int(NR2_iden[i]) << "%" << endl;
    }
    else out2b << " *" << endl;
  }
  out2b.close();

  cout << "writing clustering information" << endl;
  // write clstr information
//  I mask following 3 lines, because it crash when clusters NR
//  I thought maybe there is not a big block memory now, so
//  move the new statement to the begining of program, but because I
//  don't know the NR90_no, I just use DB_no instead
//  int *Clstr_no, *(*Clstr_list);
//  if ((Clstr_no   = new int[NR90_no]) == NULL) bomb_error("Memory");
//  if ((Clstr_list = new int*[NR90_no]) == NULL) bomb_error("Memory");

  for (i=0; i<NR90_no; i++) { Clstr_no[i]=0; Clstr_no_db1[i]=0; }
  for (i=0; i<NR_no; i++) {
    j1 = NR_clstr_no[i];
    if ( j1 < 0 ) j1 =-j1-1;
    Clstr_no[j1]++;
    Clstr_no_db1[j1]++;
  }
  for (i=0; i<NR2_no; i++) {
    j1 = NR2_clstr_no[i];
    if ( j1 >= 0) continue; // skip nr seq from db2
    if ( j1 < 0 ) j1 =-j1-1;
    Clstr_no[j1]++;
  }

  for (i=0; i<NR90_no; i++) {
    if((Clstr_list[i] = new int[ Clstr_no[i] ]) == NULL) bomb_error("Memory");
    Clstr_no[i]=0;
  }

  for (i=0; i<NR_no; i++) {
    j1 = NR_clstr_no[i];
    if ( j1 < 0 ) j1 =-j1-1;
    Clstr_list[j1][ Clstr_no[j1]++ ] = i;
  }
  for (i=0; i<NR2_no; i++) {
    j1 = NR2_clstr_no[i];
    if ( j1 >= 0) continue; // skip nr seq from db2
    if ( j1 < 0 ) j1 =-j1-1;
    Clstr_list[j1][ Clstr_no[j1]++ ] = i;
  }

  char c11;
  for (i=0; i<NR90_no; i++) {
    out2 << ">Cluster " << i << endl;
    // from db1
    for (k=0; k<Clstr_no_db1[i]; k++) {
      j = Clstr_list[i][k];
      c11 = (NR_flag[j] & IS_MINUS_STRAND) ? '-' : '+';
      out2 << k << "\t" << NR_len[j] << "nt, "<< NR_seq[j] << "...";
      if ( NR_iden[j]>0 ) 
        out2 << " at " << c11 << "/" << int(NR_iden[j]) << "%" << endl;
      else                out2 << " *" << endl;
    }

    // from db2
    for (k=Clstr_no_db1[i]; k<Clstr_no[i]; k++) {
      j = Clstr_list[i][k];
      c11 = (NR2_flag[j] & IS_MINUS_STRAND) ? '-' : '+';
      out2 << k << "\t" << NR2_len[j] << "nt, "<< NR2_seq[j] << "...";

      if ( NR2_iden[j]>0 ) {
        out2 << " at ";
        if (option_p)
          out2 << NR2_info[j][0] << ":" << NR2_info[j][1] << ":"
               << NR2_info[j][2] << ":" << NR2_info[j][3] << "/";
        out2 << c11 << "/" << int(NR2_iden[j]) << "%" << endl;
      }
      else                  out2 << " *" << endl;
    }
  }
  out2.close();
  cout << "program completed !" << endl << endl;

  times(&CPU_end);
  show_cpu_time(CPU_begin, CPU_end);

  remove_tmp_files(SEG_no, db_swap, option_B, db_bin_swap);
  remove_tmp_files_db2(option_B, db2_bin_swap);
  return 0;
} // END int main

///////////////////////FUNCTION of common tools////////////////////////////


int check_this_2d(int len, char *seqi, int &has_aas,
               int NAA, int& aan_no, int *aan_list, INTs *aan_list_no,
               INTs *look_and_count, 
               int &hit_no, int libb, int libe, int &iden_no,
               double aa1_cutoff, double aas_cutoff, double aan_cutoff,
               char this_flag, char *NR_flag, int len_upper_bound,
               int &lens, int len_lower_bound, int *check_this_info,
               double *check_this_infod) {

  static int  taap[MAX_UAA*MAX_UAA*MAX_UAA*MAX_UAA];
  static INTs aap_list[MAX_SEQ];
  static INTs aap_begin[MAX_UAA*MAX_UAA*MAX_UAA*MAX_UAA];

  int i, j, k, i1, j1, k1, i0, j0, k0, c22, sk, mm;
  int len_eff, aln_cover_flag, min_aln_lenS, min_aln_lenL;
  int required_aa1, required_aas, required_aan;

  len_eff = len;
  aln_cover_flag = 0;
  if ((check_this_infod[13] > 0.0) || (check_this_info[15]>0) ) { // has alignment coverage control
    aln_cover_flag = 1;
    min_aln_lenS = (int) (double(len) * check_this_infod[13]);
    if ( len-check_this_info[14] > min_aln_lenS)
      min_aln_lenS = len-check_this_info[14];
    if ( check_this_info[15] > min_aln_lenS)
      min_aln_lenS = check_this_info[15];
  }
  if (check_this_info[10] == 0) len_eff = min_aln_lenS; //option_G==0
  calc_required_aaxN(required_aa1, required_aas, required_aan,
                     aa1_cutoff,   aas_cutoff,   aan_cutoff, len_eff, NAA, 4);

  // check_aan_list 
  aan_no = len - NAA + 1;
  for (j=0; j<aan_no; j++) {
    aan_list[j] = 0;
    for (k=0, k1=NAA-1; k<NAA; k++, k1--) aan_list[j] += seqi[j+k] * NAAN_array[k1];
  }

  // for the short word containing 'N', mask it to '-1'
  for (j=0; j<len; j++)
    if ( seqi[j] == 4 ) {                      // here N is 4
      i0 = (j-NAA+1 > 0)      ? j-NAA+1 : 0;
      i1 = (j+NAA < aan_no)   ? j+NAA   : aan_no;
      for (i=i0; i< i1; i++) aan_list[i]=-1;
    }

  quick_sort(aan_list,0,aan_no-1);
  for(j=0; j<aan_no; j++) aan_list_no[j]=1;
  for(j=aan_no-1; j; j--) {
    if (aan_list[j] == aan_list[j-1]) {
      aan_list_no[j-1] += aan_list_no[j];
      aan_list_no[j]=0;
    }
  }
  // END check_aan_list

  // if minimal alignment length > len, return
  // I can not return earlier, because I need to calc the aan_list etc
  if (check_this_info[15]>len) return 0; // return flag=0

  // lookup_aan
  for (j=libe; j>=libb; j--) look_and_count[j]=0;
  word_table.count_word_no2(aan_no, aan_list, aan_list_no, look_and_count);


  // contained_in_old_lib()
  int band_left, band_right, best_score, band_width1, best_sum, len2, alnln, len_eff1;
  int tiden_no;
  int talign_info[5];
  int len1 = len - 4 + 1;
  char *seqj;
  int flag = 0;      // compare to old lib
  has_aas = 0;
  for (j=libb; j<=libe; j++) {
    if ( look_and_count[j] < required_aan ) continue;
    len2 = NR_len[NR90_idx[j]];
    if (len2 > len_upper_bound ) continue;
    if (len2 < len_lower_bound ) continue;
    seqj = NR_seq[NR90_idx[j]];

    if (aln_cover_flag) {
      min_aln_lenL = (int) (double(len2) * check_this_infod[11]);
      if ( len2-check_this_info[12] > min_aln_lenL)
        min_aln_lenL = len2-check_this_info[12];
      if ( check_this_info[15] > min_aln_lenL)
        min_aln_lenL = check_this_info[15];
    }
    
    if ( has_aas == 0 )  { // calculate AAP array
      for (sk=0; sk<NAA4; sk++) taap[sk] = 0;
      for (j1=0; j1<len1; j1++) {
        if ((seqi[j1]==4) || (seqi[j1+1]==4) || (seqi[j1+2]==4) || (seqi[j1+3]==4)) continue; //skip N
        c22 = seqi[j1]*NAA3 + seqi[j1+1]*NAA2 + seqi[j1+2]*NAA1 + seqi[j1+3];
        taap[c22]++;
      }
      for (sk=0,mm=0; sk<NAA4; sk++) {
        aap_begin[sk] = mm; mm+=taap[sk]; taap[sk] = 0;
      }
      for (j1=0; j1<len1; j1++) {
        if ((seqi[j1]==4) || (seqi[j1+1]==4) || (seqi[j1+2]==4) || (seqi[j1+3]==4)) continue; //skip N
        c22 = seqi[j1]*NAA3 + seqi[j1+1]*NAA2 + seqi[j1+2]*NAA1 + seqi[j1+3];
        aap_list[aap_begin[c22]+taap[c22]++] =j1;
      }
      has_aas = 1;
    }

    band_width1 = (option_b < len+len2-2 ) ? option_b : len+len2-2;
    diag_test_aapn_est(NAA1, seqj, len, len2, taap, aap_begin, 
                   aap_list, best_sum,
                   band_width1, band_left, band_right, required_aa1);
    if ( best_sum < required_aas ) continue;
    
    if (check_this_info[0] || aln_cover_flag) //return overlap region
      local_band_align2(seqi, seqj, len, len2, mat,
                        best_score, tiden_no, band_left, band_right,
                        talign_info[1],talign_info[2],
                        talign_info[3],talign_info[4], alnln);
    else
      local_band_align(seqi, seqj, len, len2, mat,
                             best_score, tiden_no, band_left, band_right);
    if ( tiden_no < required_aa1 ) continue;
    lens = (len <= len2) ? len : len2;
    len_eff1 = (check_this_info[10] == 0) ? alnln : lens;
    tiden_no = tiden_no * 100 / len_eff1;
    if (tiden_no < NR_clstr100) continue;
    if (tiden_no <= iden_no) continue; // existing iden_no
    if (aln_cover_flag) {
      if ( talign_info[4]-talign_info[3]+1 < min_aln_lenL) continue;
      if ( talign_info[2]-talign_info[1]+1 < min_aln_lenS) continue;
    }
    flag = 1; iden_no = tiden_no; hit_no = j;
    check_this_info[1] = talign_info[1];
    check_this_info[2] = talign_info[2];
    check_this_info[3] = talign_info[3];
    check_this_info[4] = talign_info[4];
    if (! check_this_info[20]) break; // not option_g
  }
  // END contained_in_old_lib()

  // comparison complimentary strand
  if (flag == 0 && option_r ) {
    for (j0=0; j0<aan_no; j0++) {
      j = aan_list[j0];
      if ( j<0 ) aan_list_comp[j0] = j;
      else       aan_list_comp[j0] = Comp_AAN_idx[j];
    }

    // lookup_aan
    for (j=libe; j>=libb; j--) look_and_count[j]=0;
    word_table.count_word_no2(aan_no, aan_list_comp, aan_list_no,
                              look_and_count);
    make_comp_iseq(len, seqi_comp, seqi);

    // reset has_aas, it use same array taap, aap_begin, aap_list
    // to store comp strand
    has_aas = 0;
    // compare to old lib
    for (j=libb; j<=libe; j++) {
      if ( look_and_count[j] < required_aan ) continue;
      len2 = NR_len[NR90_idx[j]];
      if (len2 > len_upper_bound ) continue;
      if (len2 < len_lower_bound ) continue;
      seqj = NR_seq[NR90_idx[j]];

      if (aln_cover_flag) {
        min_aln_lenL = (int) (double(len2) * check_this_infod[11]);
        if ( len2-check_this_info[12] > min_aln_lenL)
          min_aln_lenL = len2-check_this_info[12];
        if ( check_this_info[15] > min_aln_lenL)
          min_aln_lenL = check_this_info[15];
      }

      if ( has_aas == 0 )  { // calculate AAP array
        for (sk=0; sk<NAA4; sk++) taap[sk] = 0;
        for (j1=0; j1<len1; j1++) {
          if ((seqi_comp[j1]==4) || (seqi_comp[j1+1]==4) || (seqi_comp[j1+2]==4) || (seqi_comp[j1+3]==4)) continue; //skip N
          c22 = seqi_comp[j1]*NAA3 + seqi_comp[j1+1]*NAA2 + seqi_comp[j1+2]*NAA1 + seqi_comp[j1+3];
          taap[c22]++;
        }
        for (sk=0,mm=0; sk<NAA4; sk++) {
          aap_begin[sk] = mm; mm+=taap[sk]; taap[sk] = 0;
        }
        for (j1=0; j1<len1; j1++) {
          if ((seqi_comp[j1]==4) || (seqi_comp[j1+1]==4) || (seqi_comp[j1+2]==4) || (seqi_comp[j1+3]==4)) continue; //skip N
          c22 = seqi_comp[j1]*NAA3 + seqi_comp[j1+1]*NAA2 + seqi_comp[j1+2]*NAA1 + seqi_comp[j1+3];
          aap_list[aap_begin[c22]+taap[c22]++] =j1;
        }
        has_aas = 1;
      }

      band_width1 = (option_b < len+len2-2 ) ? option_b : len+len2-2;
      diag_test_aapn_est(NAA1, seqj, len, len2, taap, aap_begin,
                     aap_list, best_sum,
                     band_width1, band_left, band_right, required_aa1);
      if ( best_sum < required_aas ) continue;

      if (check_this_info[0] || aln_cover_flag) { //return overlap region
        local_band_align2(seqi_comp, seqj, len, len2, mat,
                          best_score, tiden_no, band_left, band_right,
                          talign_info[1],talign_info[2],
                          talign_info[3],talign_info[4], alnln);
        talign_info[1] = len - talign_info[1] - 1;
        talign_info[2] = len - talign_info[2] - 1;
      }
      else
        local_band_align(seqi_comp, seqj, len, len2, mat,
                         best_score, tiden_no, band_left, band_right);
      if ( tiden_no < required_aa1 ) continue;
      lens = (len <= len2) ? len : len2;
      len_eff1 = (check_this_info[10] == 0) ? alnln : lens;
      tiden_no = tiden_no * 100 / len_eff1;
      if (tiden_no < NR_clstr100) continue;
      if (tiden_no <= iden_no) continue; // existing iden_no
      if (aln_cover_flag) {
        if ( talign_info[4]-talign_info[3]+1 < min_aln_lenL) continue;
        if ( talign_info[1]-talign_info[2]+1 < min_aln_lenS) continue; //reverse
      }
      flag = -1; iden_no = tiden_no; hit_no = j;
      check_this_info[1] = talign_info[1];
      check_this_info[2] = talign_info[2];
      check_this_info[3] = talign_info[3];
      check_this_info[4] = talign_info[4];
      if (! check_this_info[20]) break; // not option_g
    }
  }
  // END if (flag == 0 && option_r )
  return flag;
} // END check_this_2d


int calc_ann_list(int len, char *seqi,
                  int NAA, int& aan_no, int *aan_list, INTs *aan_list_no) {
  int i, j, k, i1, j1, k1, i0, j0, k0, c22, sk, mm;

  // check_aan_list 
  aan_no = len - NAA + 1;
  for (j=0; j<aan_no; j++) {
    aan_list[j] = 0;
    for (k=0, k1=NAA-1; k<NAA; k++, k1--) aan_list[j] += seqi[j+k] * NAAN_array[k1];
  }

  // for the short word containing 'N', mask it to '-1'
  for (j=0; j<len; j++)
    if ( seqi[j] == 4 ) {                      // here N is 4
      i0 = (j-NAA+1 > 0)      ? j-NAA+1 : 0;
      i1 = (j+NAA < aan_no)   ? j+NAA   : aan_no;
      for (i=i0; i< i1; i++) aan_list[i]=-1;
    }

  quick_sort(aan_list,0,aan_no-1);
  for(j=0; j<aan_no; j++) aan_list_no[j]=1;
  for(j=aan_no-1; j; j--) {
    if (aan_list[j] == aan_list[j-1]) {
      aan_list_no[j-1] += aan_list_no[j];
      aan_list_no[j]=0;
    }
  }
  // END check_aan_list

  return OK_FUNC;
} // END calc_ann_list


void make_comp_iseq(int len, char *iseq_comp, char *iseq) {
  int i, j, k;
  int c[5] = {3,2,1,0,4};
  for (i=0; i<len; i++) iseq_comp[i] = c[ iseq[len-i-1] ];
} // make_comp_iseq


/////////////////////////// END ALL ////////////////////////
