/*
===============================================================================
Authors
===============================================================================

Per-Gunnar Martinsson
  Dept. of Applied Mathematics, 
  University of Colorado at Boulder, 
  526 UCB, Boulder, CO 80309-0526, USA

Gregorio Quintana-Orti
  Depto. de Ingenieria y Ciencia de Computadores, 
  Universitat Jaume I, 
  12.071 Castellon, Spain

Nathan Heavner
  Dept. of Applied Mathematics, 
  University of Colorado at Boulder, 
  526 UCB, Boulder, CO 80309-0526, USA

Robert van de Geijn
  Dept. of Computer Science and Institute for Computational Engineering and 
  Sciences, 
  The University of Texas at Austin
  Austin, TX.

===============================================================================
Copyright
===============================================================================

Copyright (C) 2016, 
  Universitat Jaume I,
  University of Colorado at Boulder,
  The University of Texas at Austin.

===============================================================================
Disclaimer
===============================================================================

This code is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY EXPRESSED OR IMPLIED.

*/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "FLAME.h"
#include "FLA_HQRRP_UT_blk_var2.h"


/*
// Matrices with dimensions smaller than THRESHOLD_FOR_DGEQPF are processed 
// with LAPACK's routine dgeqpf.
// Matrices with dimensions between THRESHOLD_FOR_DGEQPF and 
// THRESHOLD_FOR_DGEQP3 are processed with LAPACK's routine dgeqp3.
// Matrices with dimensions larger than THRESHOLD_FOR_DGEQP3 are processed 
// with the new HQRRP code.
#define THRESHOLD_FOR_DGEQPF   250
#define THRESHOLD_FOR_DGEQP3  1000
*/


// ============================================================================
// Compilation declarations.

#undef CHECK_DOWNDATING_OF_Y
#undef PROFILE


// ============================================================================
// Declaration of local prototypes.

static int MyFLA_Normal_random_matrix( FLA_Obj A ); 

static double MyFLA_Normal_random_number( double mu, double sigma );

static FLA_Error MyFLA_Downdate_Y( FLA_Obj U11, FLA_Obj U21, FLA_Obj A12,
                     FLA_Obj T, 
                     FLA_Obj Y2, FLA_Obj G1, FLA_Obj G2 );

static FLA_Error MyFLA_Apply_Q_UT_lhfc_blk_var2( FLA_Obj U11, FLA_Obj U21,
                     FLA_Obj T, 
                     FLA_Obj B1, FLA_Obj B2 );

static FLA_Error MyFLA_Apply_Q_UT_rnfc_blk_var2( FLA_Obj U11, FLA_Obj U21,
                     FLA_Obj T, 
                     FLA_Obj B1, FLA_Obj B2 );

static FLA_Error MyFLA_QRPmod_UT_unb_var2( int pivoting, int num_stages, 
                     FLA_Obj A, FLA_Obj p, FLA_Obj t, 
                     int pivot_B, FLA_Obj B,
                     int pivot_C, FLA_Obj C,
                     int build_T, FLA_Obj T );

static FLA_Error MyFLA_Apply_H2_UT_l_opd_var2( 
                     int m_u2_A2,
                     int n_a1t,
                     double * tau,
                     double * u2, int inc_u2,
                     double * a1t, int inc_a1t,
                     double * A2, int ldim_A2,
                     double * workspace );

static FLA_Error MyFLA_QRP_compute_norms( FLA_Obj A, FLA_Obj d, FLA_Obj e );

static FLA_Error NoFLA_QRP_pivot_B_C( int j_max_col,
                     int m_G, double * buff_G, int ldim_G, 
                     int pivot_B, int m_B, double * buff_B, int ldim_B, 
                     int pivot_C, int m_C, double * buff_C, int ldim_C, 
                     int * buff_p,
                     double * buff_d, double * buff_e );

static FLA_Error NoFLA_QRP_downdate_partial_norms( int m_A, int n_A,
                     double * buff_d,  int st_d,
                     double * buff_e,  int st_e,
                     double * buff_wt, int st_wt,
                     double * buff_A,  int ldim_A );

static int NoFLA_my_idamax( int n, double * buff_data );


// ============================================================================
int FLA_HQRRP_UT_blk_var2( FLA_Obj A, FLA_Obj p, FLA_Obj s, 
        int nb_alg, int pp, int panel_pivoting ) {
//
// HQRRP: It computes the Householder QR with Randomized Pivoting of matrix A.
// This routine is almost compatible with LAPACK's dgeqp3.
// The main difference is that this routine does not manage fixed columns.
//
// Main features:
//   * BLAS-3 based.
//   * Norm downdating method by Drmac.
//   * Downdating for computing Y.
//   * Use of libflame.
//   * UT transformations are used.
//
// Arguments:
// ----------
// A:              Matrix to be factorized.
// p:              output vector with the pivots.
// s:              Output vector with the tau values of the Householder factors.
// nb_alg:         Block size. Usual values for nb_alg are 32, 64, etc.
// pp:             Oversampling size. Usual values for pp are 5, 10, etc.
// panel_pivoting: If panel_pivoting==1, QR with pivoting is applied to 
//                 factorize the panels of matrix A. Otherwise, QR without 
//                 pivoting is used. Usual value for panel_pivoting is 1.
//
  // Declaration of variables.
  FLA_Obj ATL, ATR,    A00, A01, A02, 
          ABL, ABR,    A10, A11, A12,
                       A20, A21, A22;
  FLA_Obj pT,          p0,
          pB,          p1,
                       p2;
  FLA_Obj sT,          s0,
          sB,          s1,
                       s2;
  FLA_Obj GL,  GR,     G0,   G1,  G2;
  FLA_Obj YL,  YR,     Y0,   Y1,  Y2;
  FLA_Obj VL,  VR,     V0,   V1,  V2;
  FLA_Obj WL,  WR,     W0,   W1,  W2;
  FLA_Obj AR, AB1, G, Y, V, W, T1_T, None;
  int     bRow, m_A, n_A, mn_A, dtype_A, last_iter;
#ifdef PROFILE
  double  t1, t2, tt_qrp1, tt_qrp2, tt_updt, tt_down;
#endif

  // Executable Statements.
  //// printf( "%% FLA_HQRRP_UT_blk_var2.\n" );

  // Some initializations.
  m_A     = FLA_Obj_length( A );
  n_A     = FLA_Obj_width ( A );
  mn_A    = min( m_A, n_A );
  dtype_A = FLA_Obj_datatype( A );

  // Quick return.
  if( mn_A == 0 ) {
    return 0;
  }

  // Initialize the seed for the generator of random numbers.
  srand( 12 );

#ifdef PROFILE
  tt_qrp1 = 0.0;
  tt_qrp2 = 0.0;
  tt_updt = 0.0;
  tt_down = 0.0;
#endif

  // Create and initialize auxiliary objects.
  FLA_Obj_create( dtype_A, nb_alg + pp, m_A,    0, 0, & G );
  FLA_Obj_create( dtype_A, nb_alg + pp, n_A,    0, 0, & Y );
  FLA_Obj_create( dtype_A, nb_alg + pp, n_A,    0, 0, & V );
  FLA_Obj_create( dtype_A, nb_alg,      n_A,    0, 0, & W );

  // Initialize matrices G and Y.
  MyFLA_Normal_random_matrix( G );
  FLA_Gemm( FLA_NO_TRANSPOSE, FLA_NO_TRANSPOSE, 
            FLA_ONE, G, A, FLA_ZERO, Y );

  // Initial Partitioning.
  FLA_Part_2x2( A,    & ATL, & ATR,
                      & ABL, & ABR,    0, 0, FLA_TL );
  FLA_Part_2x1( p,    & pT, 
                      & pB,            0, FLA_TOP );
  FLA_Part_2x1( s,    & sT,
                      & sB,            0, FLA_TOP );
  FLA_Part_1x2( G,    & GL,  & GR,     0, FLA_LEFT );
  FLA_Part_1x2( Y,    & YL,  & YR,     0, FLA_LEFT );
  FLA_Part_1x2( V,    & VL,  & VR,     0, FLA_LEFT );
  FLA_Part_1x2( W,    & WL,  & WR,     0, FLA_LEFT );

  // Main Loop.
  while( ( FLA_Obj_length( ATL ) < FLA_Obj_length( A ) )&&
         ( FLA_Obj_width ( ATL ) < FLA_Obj_width ( A ) ) ) {
    bRow = min( FLA_Obj_min_dim( ABR ), nb_alg );

    // Iteration Initial Partitioning.
    FLA_Repart_2x2_to_3x3( ATL, /**/ ATR,       &A00, /**/ &A01, &A02,
                        /* ************* */   /* ******************** */
                                                &A10, /**/ &A11, &A12,
                           ABL, /**/ ABR,       &A20, /**/ &A21, &A22,
                           bRow, bRow, FLA_BR );
    FLA_Repart_2x1_to_3x1( pT,                  &p0, 
                        /* ** */              /* *** */
                                                &p1, 
                           pB,                  &p2,    bRow, FLA_BOTTOM );
    FLA_Repart_2x1_to_3x1( sT,                  &s0, 
                        /* ** */              /* ** */
                                                &s1, 
                           sB,                  &s2,        bRow, FLA_BOTTOM );
    FLA_Repart_1x2_to_1x3( GL,   /**/ GR,       &G0,  /**/ &G1,  &G2,
                           bRow, FLA_RIGHT );
    FLA_Repart_1x2_to_1x3( YL,   /**/ YR,       &Y0,  /**/ &Y1,  &Y2,
                           bRow, FLA_RIGHT );
    FLA_Repart_1x2_to_1x3( VL,   /**/ VR,       &V0,  /**/ &V1,  &V2,
                           bRow, FLA_RIGHT );
    FLA_Repart_1x2_to_1x3( WL,   /**/ WR,       &W0,  /**/ &W1,  &W2,
                           bRow, FLA_RIGHT );

    // ------------------------------------------------------------------------

    //// printf( "Iter: %ld \n", FLA_Obj_length( ATL ) );

    // Check whether it is the last iteration.
    last_iter = ( FLA_Obj_min_dim( A22 ) <= 0 ? 1 : 0 );

#ifdef CHECK_DOWNDATING_OF_Y
    // Check downdating of matrix Y: Compare downdated matrix Y with 
    // matrix Y computed from scratch.
    int      m_cyr, n_cyr, ldim_cyr, ldim_YR, ii, jj; 
    double   * buff_cyr, * buff_YR, aux, sum;
    FLA_Obj  CYR;

    FLA_Obj_create_conf_to( FLA_NO_TRANSPOSE, YR, & CYR );

    FLA_Gemm( FLA_NO_TRANSPOSE, FLA_NO_TRANSPOSE, 
              FLA_ONE, GR, ABR, FLA_ZERO, CYR ); 

    //// FLA_Obj_show( " CYR = [ ", CYR, "%le", " ];" );
    //// FLA_Obj_show( " Y  = [ ", Y, "%le", " ];" );

    m_cyr    = FLA_Obj_length( CYR );
    n_cyr    = FLA_Obj_width( CYR );
    buff_cyr = ( double * ) FLA_Obj_buffer_at_view( CYR );
    ldim_cyr = FLA_Obj_col_stride( CYR );
    buff_YR   = ( double * ) FLA_Obj_buffer_at_view( YR );
    ldim_YR  = FLA_Obj_col_stride( YR );
    sum = 0.0;
    for( jj = 0; jj < n_cyr; jj++ ) {
      for( ii = 0; ii < m_cyr; ii++ ) {
        aux = buff_YR [ ii + jj * ldim_YR ] -
              buff_cyr[ ii + jj * ldim_cyr ];
        sum += aux * aux;
      }
    }
    sum = sqrt( sum );
    printf( "%%  diff between Y and downdated Y: %le\n", sum );

    FLA_Obj_free( & CYR );
#endif

    if( last_iter == 0 ) {
      // Compute QRP of YR, and apply permutations to matrix AR.
      // A copy of YR is made into VR, and permutations are applied to YR.
#ifdef PROFILE
      t1 = FLA_Clock();
#endif
      FLA_Merge_2x1( ATR,
                     ABR,   & AR );
      FLA_Copy( YR, VR );
      MyFLA_QRPmod_UT_unb_var2( 1, bRow, VR, pB, sB, 1, AR, 1, YR, 0, None );

#ifdef PROFILE
      t2 = FLA_Clock();
      tt_qrp1 += ( t2 - t1 );
#endif
    }

    //
    // Compute QRP of panel AB1 = [ A11; A21 ].
    // Apply same permutations to A01 and Y1, and build T1_T.
    //
#ifdef PROFILE
    t1 = FLA_Clock();
#endif

    FLA_Part_2x1( W1,   & T1_T,
                        & None,    bRow, FLA_TOP );
    FLA_Merge_2x1( A11,
                   A21,   & AB1 );
    MyFLA_QRPmod_UT_unb_var2( panel_pivoting, -1, AB1, p1, s1, 
        1, A01, 1, Y1, 1, T1_T );

#ifdef PROFILE
    t2 = FLA_Clock();
    tt_qrp2 += ( t2 - t1 );
#endif

    //
    // Update the rest of the matrix.
    //
    if ( FLA_Obj_width( A12 ) > 0 ) {
      // Apply the Householder transforms associated with AB1 = [ A11; A21 ] 
      // and T1_T to [ A12; A22 ]:
      //   / A12 \ := QB1' / A12 \
      //   \ A22 /         \ A22 /
      // where QB1 is formed from AB1 and T1_T.
#ifdef PROFILE
      t1 = FLA_Clock();
#endif

      MyFLA_Apply_Q_UT_lhfc_blk_var2( A11, A21, T1_T, A12, A22 );

#ifdef PROFILE
      t2 = FLA_Clock();
      tt_updt += ( t2 - t1 );
#endif
    }

    //
    // Downdate matrix Y.
    //
    if ( FLA_Obj_width( Y2 ) > 0 ) {
#ifdef PROFILE
      t1 = FLA_Clock();
#endif

      MyFLA_Downdate_Y( A11, A21, A12, T1_T, Y2, G1, G2 );

#ifdef PROFILE
      t2 = FLA_Clock();
      tt_down += ( t2 - t1 );
#endif
    }

    // ------------------------------------------------------------------------
    // Iteration Final Repartitioning.

    FLA_Cont_with_3x3_to_2x2( &ATL, /**/ &ATR,       A00, A01, /**/ A02,
                                                     A10, A11, /**/ A12,
                            /* ************** */  /* ****************** */
                              &ABL, /**/ &ABR,       A20, A21, /**/ A22,
                              FLA_TL );
    FLA_Cont_with_3x1_to_2x1( &pT,                   p0, 
                                                     p1, 
                            /* ** */              /* *** */
                              &pB,                   p2,     FLA_TOP );
    FLA_Cont_with_3x1_to_2x1( &sT,                   s0, 
                                                     s1, 
                            /* ** */              /* ** */
                              &sB,                   s2,     FLA_TOP );
    FLA_Cont_with_1x3_to_1x2( &GL,   /**/ &GR,       G0,  G1,  /**/ G2,
                              FLA_LEFT );
    FLA_Cont_with_1x3_to_1x2( &YL,   /**/ &YR,       Y0,  Y1,  /**/ Y2,
                              FLA_LEFT );
    FLA_Cont_with_1x3_to_1x2( &VL,   /**/ &VR,       V0,  V1,  /**/ V2,
                              FLA_LEFT );
    FLA_Cont_with_1x3_to_1x2( &WL,   /**/ &WR,       W0,  W1,  /**/ W2,
                              FLA_LEFT );
  }

  // Remove auxiliary objects.
  FLA_Obj_free( & G );
  FLA_Obj_free( & Y );
  FLA_Obj_free( & V );
  FLA_Obj_free( & W );

#ifdef PROFILE
  printf( "%% tt_qrp_1:       %le\n", tt_qrp1 );
  printf( "%% tt_qrp_2:       %le\n", tt_qrp2 );
  printf( "%% tt_updt_rest:   %le\n", tt_updt );
  printf( "%% tt_downdating:  %le\n", tt_down );
  printf( "%% total_time:     %le\n", tt_qrp1 + tt_qrp2 + tt_updt + tt_down);
#endif

  return 0;
}


// ============================================================================
static int MyFLA_Normal_random_matrix( FLA_Obj A ) {
//
// It sets a random matrix with normal distribution.
//

  switch( FLA_Obj_datatype( A ) ) {

  case FLA_DOUBLE:
    {
      double  * buff_A;
      int     m_A, n_A, ldim_A, i, j;

      // Some initializations.
      m_A     = FLA_Obj_length( A );
      n_A     = FLA_Obj_width ( A );
      buff_A  = ( double * ) FLA_Obj_buffer_at_view( A );
      ldim_A  = FLA_Obj_col_stride( A );

      // Main loop.
      for ( j = 0; j < n_A; j++ ) {
        for ( i = 0; i < m_A; i++ ) {
          buff_A[ i + j * ldim_A ] = MyFLA_Normal_random_number( 0.0, 1.0 );
        }
      }
    }
    break;

  default:
    fprintf( stderr, "+++ ERROR in MyFLA_Normal_random_matrix:  " );
    fprintf( stderr, "Datatype not implemented: %d\n", FLA_Obj_datatype( A ) );
  }

  return FLA_SUCCESS;
}

/* ========================================================================= */
static double MyFLA_Normal_random_number( double mu, double sigma ) {
  static int     alternate_calls = 0;
  static double  b1, b2;
  double         c1, c2, a, factor;

  // Quick return.
  if( alternate_calls == 1 ) {
    alternate_calls = ! alternate_calls;
    return( mu + sigma * b2 );
  }
  // Main loop.
  do {
    c1 = -1.0 + 2.0 * ( (double) rand() / RAND_MAX );
    c2 = -1.0 + 2.0 * ( (double) rand() / RAND_MAX );
    a = c1 * c1 + c2 * c2;
  } while ( ( a == 0 )||( a >= 1 ) );
  factor = sqrt( ( -2 * log( a ) ) / a );
  b1 = c1 * factor;
  b2 = c2 * factor;
  alternate_calls = ! alternate_calls;
  return( mu + sigma * b1 );
}

// ============================================================================
static FLA_Error MyFLA_Downdate_Y( FLA_Obj U11, FLA_Obj U21, FLA_Obj A12,
                     FLA_Obj T, 
                     FLA_Obj Y2, FLA_Obj G1, FLA_Obj G2 ) {
//
// It downdates matrix Y, and updates matrix G.
// Only Y2 of Y is updated.
// Only G1 and G2 of G are updated.
//
// Y2 = Y2 - ( G1 - ( G1*U11 + G2*U21 ) * inv( T11 ) * U11' ) * R12.
//
  FLA_Obj  B;

  // Create auxiliary object B.
  FLA_Obj_create_conf_to( FLA_NO_TRANSPOSE, G1, & B );

  // B = G1.
  FLA_Copy( G1, B );

  // B = B * U11.
  FLA_Trmm( FLA_RIGHT, FLA_LOWER_TRIANGULAR,
            FLA_NO_TRANSPOSE, FLA_UNIT_DIAG,
            FLA_ONE, U11, B );

  // B = B + G2 * U21.
  FLA_Gemm( FLA_NO_TRANSPOSE, FLA_NO_TRANSPOSE,
            FLA_ONE, G2, U21, FLA_ONE, B );

  // B = B * inv( T11 ).
  FLA_Trsm( FLA_RIGHT, FLA_UPPER_TRIANGULAR,
            FLA_NO_TRANSPOSE, FLA_NONUNIT_DIAG,
            FLA_ONE, T, B );

  // B = - B * U11^H.
  FLA_Trmm( FLA_RIGHT, FLA_LOWER_TRIANGULAR,
            FLA_CONJ_TRANSPOSE, FLA_UNIT_DIAG,
            FLA_MINUS_ONE, U11, B );

  // B = G1 + B.
  FLA_Axpy( FLA_ONE, G1, B );

  // Y2 = Y2 - B * R12.
  FLA_Gemm( FLA_NO_TRANSPOSE, FLA_NO_TRANSPOSE,
            FLA_MINUS_ONE, B, A12, FLA_ONE, Y2 );

  //
  // GR = GR * Q
  //
  MyFLA_Apply_Q_UT_rnfc_blk_var2( U11, U21, T, G1, G2 );

  // Remove auxiliary object B.
  FLA_Obj_free( & B );

  return FLA_SUCCESS;
}

// ============================================================================
static FLA_Error MyFLA_Apply_Q_UT_lhfc_blk_var2( FLA_Obj U11, FLA_Obj U21,
                     FLA_Obj T, 
                     FLA_Obj B1, FLA_Obj B2 ) {
//
// It applies the conjugate-transpose of a unitary matrix Q to a matrix B from 
// the left:
//   B := Q' B
// where:
//   B = [ B1; B2 ]
//   U = [ U11; U21 ]
//   Q = ( I - U * inv(T)' * U' )'.
//
  FLA_Obj  W;

  // Create auxiliary object.
  FLA_Obj_create_conf_to( FLA_NO_TRANSPOSE, B1, & W );

  // W = B1;

  FLA_Copyt( FLA_NO_TRANSPOSE, B1, W );

  // U11 = trilu( U11 );
  // U21 = U21;
  // W = triu( T )' * ( U11' * B1 + U21' * B2 );

  FLA_Trmm( FLA_LEFT, FLA_LOWER_TRIANGULAR,
            FLA_CONJ_TRANSPOSE, FLA_UNIT_DIAG,
            FLA_ONE, U11, W );

  FLA_Gemm( FLA_CONJ_TRANSPOSE, FLA_NO_TRANSPOSE, 
            FLA_ONE, U21, B2, FLA_ONE, W );

  FLA_Trsm( FLA_LEFT, FLA_UPPER_TRIANGULAR,
            FLA_CONJ_TRANSPOSE, FLA_NONUNIT_DIAG,
            FLA_ONE, T, W );

  // B2 = B2 - U21 * W;
  // B1 = B1 - U11 * W;

  FLA_Gemm( FLA_NO_TRANSPOSE, FLA_NO_TRANSPOSE,
            FLA_MINUS_ONE, U21, W, FLA_ONE, B2 );

  FLA_Trmm( FLA_LEFT, FLA_LOWER_TRIANGULAR,
            FLA_NO_TRANSPOSE, FLA_UNIT_DIAG,
            FLA_MINUS_ONE, U11, W );

  FLA_Axpyt( FLA_NO_TRANSPOSE, FLA_ONE, W, B1 );

  // Remove auxiliary object.
  FLA_Obj_free( & W );

  return FLA_SUCCESS;
}

// ============================================================================
static FLA_Error MyFLA_Apply_Q_UT_rnfc_blk_var2( FLA_Obj U11, FLA_Obj U21, 
                     FLA_Obj T,
                     FLA_Obj B1, FLA_Obj B2 ) {
//
// It applies a unitary matrix Q to a matrix B from the right:
//   B = B * Q
// where:
//   B = [ B1; B2 ]
//   U = [ U11; U21 ]
//   Q = ( I - U * inv(T)' * U' ).
//
  FLA_Obj  W;

  // Create auxiliary object.
  FLA_Obj_create_conf_to( FLA_TRANSPOSE, B1, & W );

  // W = B1^T;

  FLA_Copyt( FLA_TRANSPOSE, B1, W );

  // U11 = trilu( U11 );
  // U21 = U21;
  // Let W^T be conformal to B1.
  // W^T = ( B1 * U11 + B2 * U21 ) * inv( triu(T) );
  // W   = triu(T)' * ( U11' * B1' + U21' * B2' );

  FLA_Trmm( FLA_LEFT, FLA_LOWER_TRIANGULAR,
            FLA_TRANSPOSE, FLA_UNIT_DIAG,
            FLA_ONE, U11, W );

  FLA_Gemm( FLA_TRANSPOSE, FLA_TRANSPOSE, 
            FLA_ONE, U21, B2, FLA_ONE, W );

  FLA_Trsm( FLA_LEFT, FLA_UPPER_TRIANGULAR,
            FLA_TRANSPOSE, FLA_NONUNIT_DIAG,
            FLA_ONE, T, W );

  // B2 = B2 - W^T * U21';
  // B1 = B1 - W^T * U11';
  //    = B1 - ( conj(U11) * W )^T;

  FLA_Gemm( FLA_TRANSPOSE, FLA_CONJ_TRANSPOSE,
            FLA_MINUS_ONE, W, U21, FLA_ONE, B2 );

  FLA_Trmm( FLA_LEFT, FLA_LOWER_TRIANGULAR,
            FLA_CONJ_NO_TRANSPOSE, FLA_UNIT_DIAG,
            FLA_MINUS_ONE, U11, W );

  FLA_Axpyt( FLA_TRANSPOSE, FLA_ONE, W, B1 );

  // Remove auxiliary object.
  FLA_Obj_free( & W );

  return FLA_SUCCESS;
}

/* ========================================================================= */
static FLA_Error MyFLA_QRPmod_UT_unb_var2( int pivoting, int num_stages, 
                     FLA_Obj A, FLA_Obj p, FLA_Obj t, 
                     int pivot_B, FLA_Obj B,
                     int pivot_C, FLA_Obj C,
                     int build_T, FLA_Obj T ) {
//
// It computes an unblocked QR factorization of matrix A with or without 
// pivoting. Matrices B and C are optionally pivoted, and matrix T is
// optionally built.
//
// Arguments:
// "pivoting": If pivoting==1, then QR factorization with pivoting is used.
// "numstages": It tells the number of columns that are factorized.
//   If "num_stages" is negative, the whole matrix A is factorized.
//   If "num_stages" is positive, only the first "num_stages" are factorized.
// "pivot_B": if "pivot_B" is true, matrix "B" is pivoted too.
// "pivot_C": if "pivot_C" is true, matrix "C" is pivoted too.
// "build_T": if "build_T" is true, matrix "T" is built.
//
  // Declaration of variables.
  FLA_Obj  d, e, workspace;
  int      j, m_A, n_A, mn_A, m_B, m_C, dtype_A, ldim_A, ldim_B, ldim_C, 
           ldim_T, m_A20, n_A20, m_a21, m_A22, n_A22, n_dB, idx_max_col, 
           i_one = 1;
  double   * buff_A, * buff_t, * buff_B, * buff_C, * buff_T, 
           * buff_d, * buff_e, * buff_workspace, d_one = 1.0;
  int      * buff_p;

  int idamax_();

  //// printf( "  %% MyFLA_QRPmod_UT_unb_var2. build_T: %d\n", build_T );

  // Some initializations.
  dtype_A = FLA_Obj_datatype( A );
  m_A     = FLA_Obj_length( A );
  n_A     = FLA_Obj_width ( A );
  mn_A    = min( m_A, n_A );
  ldim_A  = FLA_Obj_col_stride( A );
  buff_A  = ( double * ) FLA_Obj_buffer_at_view( A );

  buff_p  = ( int * )    FLA_Obj_buffer_at_view( p );

  buff_t  = ( double * ) FLA_Obj_buffer_at_view( t );

  if( pivot_B ) {
    m_B     = FLA_Obj_length( B );
    buff_B  = ( double * ) FLA_Obj_buffer_at_view( B );
    ldim_B  = FLA_Obj_col_stride( B );
  }

  if( pivot_C ) {
    m_C     = FLA_Obj_length( C );
    buff_C  = ( double * ) FLA_Obj_buffer_at_view( C );
    ldim_C  = FLA_Obj_col_stride( C );
  }

  if( build_T ) {
    buff_T  = ( double * ) FLA_Obj_buffer_at_view( T );
    ldim_T  = FLA_Obj_col_stride( T );
  }

  // Set the number of stages, if needed.
  if( num_stages < 0 ) {
    num_stages = mn_A;
  }

  // Create auxiliary objects.
  FLA_Obj_create( dtype_A, n_A, 1, 0, 0, & d );
  FLA_Obj_create( dtype_A, n_A, 1, 0, 0, & e );
  FLA_Obj_create( dtype_A, n_A, 1, 0, 0, & workspace );

  buff_d = ( double * ) FLA_Obj_buffer_at_view( d );
  buff_e = ( double * ) FLA_Obj_buffer_at_view( e );
  buff_workspace = ( double * ) FLA_Obj_buffer_at_view( workspace );

  if( pivoting == 1 ) {
    // Compute initial norms of A into d and e.
    MyFLA_QRP_compute_norms( A, d, e );
  }

  // Main Loop.
  //// printf( "  m x n: %d x %d  num_stages: %d \n", m_A, n_A, num_stages );
  for( j = 0; j < num_stages; j++ ) {
    //// printf( "%% ----------------------------------\n" );
    //// printf( "%% Iter: %d \n", j );
    //// printf( "%% ----------------------------------\n" );

    n_dB  = n_A - j;
    m_a21 = m_A - j - 1;
    m_A22 = m_A - j - 1;
    n_A22 = n_A - j - 1;
    m_A20 = m_A - j - 1;
    n_A20 = j;

    if( pivoting == 1 ) {
      // Obtain the index of the column with largest 2-norm.
      //// FLA_Amax( dB, max_col );
      // Function idamax of MKL can fail by returning an index equal or larger 
      // than the dimension when working with small numbers. Hence, an own 
      // implementation is used.
      //// printf( "  Chivato 901. n_dB: %d  my_idamax:   %d\n", 
      ////         n_dB, NoFLA_my_idamax( n_dB, & buff_d[ j ] ) );
      //// printf( "  Chivato 902. n_dB: %d  idx_max_col: %d\n", 
      ////         n_dB, idamax_( & n_dB, & buff_d[ j ], & i_one ) - 1 );
      //// idx_max_col = idamax_( & n_dB, & buff_d[ j ], & i_one ) - 1;
      idx_max_col = NoFLA_my_idamax( n_dB, & buff_d[ j ] );

      // Swap columns of G, pivots, and norms.
      //// FLA_QRP_pivot_B( max_col, GR, BR, pB, dB, eB );
      NoFLA_QRP_pivot_B_C( idx_max_col,
          m_A, & buff_A[ 0 + j * ldim_A ], ldim_A,
          pivot_B, m_B, & buff_B[ 0 + j * ldim_B ], ldim_B,
          pivot_C, m_C, & buff_C[ 0 + j * ldim_C ], ldim_C,
          & buff_p[ j ],
          & buff_d[ j ],
          & buff_e[ j ] );
    }

    // Compute tau1 and u21 from alpha11 and a21 such that tau1 and u21
    // determine a Householder transform H such that applying H from the
    // left to the column vector consisting of alpha11 and a21 annihilates
    // the entries in a21 (and updates alpha11).
    //// FLA_Househ2( alpha11,
    ////              a21, tau1 );
    //// FLA_Househ2_UT( FLA_LEFT,
    ////                 alpha11,
    ////                 a21, tau1 );

    FLA_Househ2_UT_l_opd( m_a21,
                          & buff_A[ j + j * ldim_A ],
                          & buff_A[ ( j+1 ) + j * ldim_A ], i_one,
                          & buff_t[ j ] );

    // / a12t \ =  H / a12t \
    // \ A22  /      \ A22  /
    //
    // where H is formed from tau1 and u21.
    //// FLA_QR_Update_Rest_blk_b2( tau1, a21, a12t,
    ////                                       A22 );
    /// FLA_Apply_H2_UT( FLA_LEFT, tau1, a21, a12t,
    ///                                       A22 );

    MyFLA_Apply_H2_UT_l_opd_var2(
        m_A22, //// m_u2_A2,
        n_A22, //// n_a1t,
        & buff_t[ j ], //// tau,
        & buff_A[ ( j+1 ) + j * ldim_A ], 1, //// u2,inc_u2,
        & buff_A[ j + ( j+1 ) * ldim_A ], ldim_A, //// a1t,inc_a1t,
        & buff_A[ ( j+1 ) + ( j+1 ) * ldim_A ], ldim_A,
        buff_workspace ); //// A2,rs_A2,cs_A2);

    // Build T.
    if( build_T ) {
      // rho11 = tau1;
      //// FLA_Copy( tau1, rho11 );
      buff_T[ j + j * ldim_T ] = buff_t[ j ];

      // t01 = a10t' + A20' * u21;
      //// FLA_Copyt_external( FLA_CONJ_TRANSPOSE, a10t, r01 );
      dcopy_( & j, & buff_A[ j + 0 * ldim_A ], & ldim_A,
                   & buff_T[ 0 + j * ldim_T ], & i_one );

      //// FLA_Gemv_external( FLA_CONJ_TRANSPOSE, FLA_ONE, A20, a21, 
      ////                    FLA_ONE, r01 );
      dgemv_( "Transpose", & m_A20, & n_A20,
          & d_one,
          & buff_A[ ( j+1 ) + 0 * ldim_A ], & ldim_A,
          & buff_A[ ( j+1 ) + j * ldim_A ], & i_one,
          & d_one,
          & buff_T[ 0 + j * ldim_T ], & i_one );
    }

    if( pivoting == 1 ) {
      // Update partial column norms.
      //// FLA_QRP_update_partial_norms( dB, eB, a12t, A22 );

      NoFLA_QRP_downdate_partial_norms( m_A22, n_A22, 
          & buff_d[ j+1 ], 1,
          & buff_e[ j+1 ], 1,
          & buff_A[ j + ( j+1 ) * ldim_A ], ldim_A,
          & buff_A[ ( j+1 ) + ( j+1 ) * ldim_A ], ldim_A );
    }
  }

  // Remove auxiliary objects.
  FLA_Obj_free( & d );
  FLA_Obj_free( & e );
  FLA_Obj_free( & workspace );

  return FLA_SUCCESS;
}

/* ========================================================================= */
static FLA_Error MyFLA_Apply_H2_UT_l_opd_var2( 
                     int m_u2_A2,
                     int n_a1t,
                     double * tau,
                     double * u2, int inc_u2,
                     double * a1t, int inc_a1t,
                     double * A2, int ldim_A2, 
                     double * workspace ) {

  double  one_p       = 1.0;
  double  minus_one_p = -1.0;
  double  rtau;
  int     inc_w1t;

  // FLA_Obj w1t;
  double * w1t;

  // if ( FLA_Obj_has_zero_dim( a1t ) ) return FLA_SUCCESS;
  if ( n_a1t == 0 ) return FLA_SUCCESS;

  // FLA_Obj_create_conf_to( FLA_NO_TRANSPOSE, a1t, &w1t );
  //// w1t = ( double* ) FLA_malloc( n_a1t * sizeof( *a1t ) );
  w1t = workspace;
  inc_w1t = 1;

  // // w1t = a1t;
  // FLA_Copy_external( a1t, w1t );
  //// bl1_dcopyv( BLIS1_NO_CONJUGATE,
  ////             n_a1t,
  ////             a1t, inc_a1t, 
  ////             w1t, inc_w1t ); 
  dcopy_( & n_a1t,
          a1t, & inc_a1t, 
          w1t, & inc_w1t ); 

  // // w1t = w1t + u2' * A2;
  // // w1t = w1t + A2^T * conj(u2);
  // FLA_Gemvc_external( FLA_TRANSPOSE, FLA_CONJUGATE, FLA_ONE, A2, u2, 
  //     FLA_ONE, w1t );
  //// bl1_dgemv( BLIS1_TRANSPOSE,
  ////            BLIS1_CONJUGATE,
  ////            m_u2_A2,
  ////            n_a1t,
  ////            one_p,
  ////            A2, rs_A2, cs_A2,
  ////            u2, inc_u2,
  ////            one_p,
  ////            w1t, inc_w1t );
  dgemv_( "Transpose",
          & m_u2_A2,
          & n_a1t,
          & one_p,
          A2, & ldim_A2,
          u2, & inc_u2,
          & one_p,
          w1t, & inc_w1t );

  // // w1t = w1t / tau;
  // FLA_Inv_scalc_external( FLA_NO_CONJUGATE, tau, w1t );
  //// bl1_dinvscalv( BLIS1_NO_CONJUGATE,
  ////                n_a1t,
  ////                tau,
  ////                w1t, inc_w1t );
  if( * tau == 0.0 ) {
    fprintf( stderr, "ERROR in MyFLA_Apply_H2_UT_l_opd_var2: Tau is zero.\n" );
  } else {
    rtau = 1.0 / ( * tau );
  }
  dscal_( & n_a1t,
          & rtau,
          w1t, & inc_w1t );

  // // a1t = - w1t + a1t;
  // FLA_Axpy_external( FLA_MINUS_ONE, w1t, a1t );
  //// bl1_daxpyv( BLIS1_NO_CONJUGATE,
  ////             n_a1t,
  ////             minus_one_p,
  ////             w1t, inc_w1t,
  ////             a1t, inc_a1t );
  daxpy_( & n_a1t,
          & minus_one_p,
          w1t, & inc_w1t,
          a1t, & inc_a1t );

  // // A2 = - u2 * w1t + A2;
  // FLA_Ger_external( FLA_MINUS_ONE, u2, w1t, A2 );
  //// bl1_dger( BLIS1_NO_CONJUGATE,
  ////           BLIS1_NO_CONJUGATE,
  ////           m_u2_A2,
  ////           n_a1t,
  ////           minus_one_p,
  ////           u2, inc_u2,
  ////           w1t, inc_w1t,
  ////           A2, rs_A2, cs_A2 );
  dger_( & m_u2_A2,
         & n_a1t,
         & minus_one_p,
         u2, & inc_u2,
         w1t, & inc_w1t,
         A2, & ldim_A2 );

  //// // FLA_Obj_free( &w1t );
  //// FLA_free( w1t );

  return FLA_SUCCESS;
}

/* ========================================================================= */
static FLA_Error MyFLA_QRP_compute_norms( FLA_Obj A, FLA_Obj d, FLA_Obj e ) {

  switch( FLA_Obj_datatype( A ) ) {

  case FLA_DOUBLE: {
    int     m_A, n_A, ld_A, st_d, st_e, j, i_one = 1;
    double  * ptr_A, * ptr_d, * ptr_e, dnrm2_();

    m_A   = FLA_Obj_length( A );
    n_A   = FLA_Obj_width( A );
    ptr_A = ( double * ) FLA_Obj_buffer_at_view( A );
    ptr_d = ( double * ) FLA_Obj_buffer_at_view( d );
    ptr_e = ( double * ) FLA_Obj_buffer_at_view( e );
    ld_A  = FLA_Obj_col_stride( A );

    st_d = 0;
    if( FLA_Obj_length( d ) == 1 ) {
      st_d = FLA_Obj_col_stride( d );
    } else if( FLA_Obj_width( d ) == 1 ) {
      st_d = 1;
    } else {
      FLA_Print_message( "MyFLA_QRP_compute_norms: Object d is not a vector", 
                         __FILE__, __LINE__ );
      FLA_Abort();
    }
    st_e = 0;
    if( FLA_Obj_length( e ) == 1 ) {
      st_e = FLA_Obj_col_stride( e );
    } else if( FLA_Obj_width( e ) == 1 ) {
      st_e = 1;
    } else {
      FLA_Print_message( "MyFLA_QRP_compute_norms: Object e is not a vector", 
                         __FILE__, __LINE__ );
      FLA_Abort();
    }

    for( j = 0; j < n_A; j++ ) {
      *ptr_d = dnrm2_( & m_A, ptr_A, & i_one );
      *ptr_e = *ptr_d;
      ptr_A += ld_A;
      ptr_d += st_d;
      ptr_e += st_e;
    }

    break;
  }
  default:
    FLA_Print_message( "MyFLA_QRP_compute_norms: datatype not yet implemented", 
                       __FILE__, __LINE__ );
    FLA_Abort();
  }

  return FLA_SUCCESS;
}

/* ========================================================================= */
static FLA_Error NoFLA_QRP_pivot_B_C( int j_max_col,
                     int m_G, double * buff_G, int ldim_G, 
                     int pivot_B, int m_B, double * buff_B, int ldim_B, 
                     int pivot_C, int m_C, double * buff_C, int ldim_C, 
                     int * buff_p,
                     double * buff_d, double * buff_e ) {
// Matrices B and C are optionally pivoted.

  // Declaration of variables.
  int     ival, i_one = 1;
  double  * ptr_g1, * ptr_g2, * ptr_b1, * ptr_b2, * ptr_c1, * ptr_c2;

  //// j_max_col = *( ( int * ) FLA_Obj_buffer_at_view( max_col ) );
  //// printf( "j_max_col: %d\n", j_max_col );
  //// FLA_Obj_show( " max_col = [ ", max_col, "%d", " ];" );

  // Swap columns of G, pivots, and norms.
  if( j_max_col != 0 ) {

    // Swap full column 0 and column "j_max_col" of G.
    ptr_g1 = & buff_G[ 0 + 0         * ldim_G ];
    ptr_g2 = & buff_G[ 0 + j_max_col * ldim_G ];
    dswap_( & m_G, ptr_g1, & i_one, ptr_g2, & i_one );

    // Swap full column 0 and column "j_max_col" of B.
    if( pivot_B ) {
      ptr_b1 = & buff_B[ 0 + 0         * ldim_B ];
      ptr_b2 = & buff_B[ 0 + j_max_col * ldim_B ];
      dswap_( & m_B, ptr_b1, & i_one, ptr_b2, & i_one );
    }

    // Swap full column 0 and column "j_max_col" of C.
    if( pivot_C ) {
      ptr_c1 = & buff_C[ 0 + 0         * ldim_C ];
      ptr_c2 = & buff_C[ 0 + j_max_col * ldim_C ];
      dswap_( & m_C, ptr_c1, & i_one, ptr_c2, & i_one );
    }

    // Swap element 0 and element "j_max_col" of pivot vector "p".
    ival = buff_p[ j_max_col ];
    buff_p[ j_max_col ] = buff_p[ 0 ];
    buff_p[ 0 ] = ival;

    // Copy norms of column 0 to column "j_max_col".
    buff_d[ j_max_col ] = buff_d[ 0 ];
    buff_e[ j_max_col ] = buff_e[ 0 ];
  }

  return FLA_SUCCESS;
}

/* ========================================================================= */
static FLA_Error NoFLA_QRP_downdate_partial_norms( int m_A, int n_A,
                     double * buff_d,  int st_d,
                     double * buff_e,  int st_e,
                     double * buff_wt, int st_wt,
                     double * buff_A,  int ldim_A ) {
  //// printf( "  NoFLA_QRP_downdate_partial_norms by Drmac. \n" );

  int     j, i_one = 1;
  double  * ptr_d, * ptr_e, * ptr_wt, * ptr_A;
  double  temp, temp2, temp5, tol3z;
  double  dnrm2_();

  /*
*
*           Update partial column norms
*
          DO 30 J = I + 1, N
             IF( WORK( J ).NE.ZERO ) THEN
*
*                 NOTE: The following 4 lines follow from the analysis in
*                 Lapack Working Note 176.
*                 
                TEMP = ABS( A( I, J ) ) / WORK( J )
                TEMP = MAX( ZERO, ( ONE+TEMP )*( ONE-TEMP ) )
                TEMP2 = TEMP*( WORK( J ) / WORK( N+J ) )**2
                IF( TEMP2 .LE. TOL3Z ) THEN 
                   IF( M-I.GT.0 ) THEN
                      WORK( J ) = DNRM2( M-I, A( I+1, J ), 1 )
                      WORK( N+J ) = WORK( J )
                   ELSE
                      WORK( J ) = ZERO
                      WORK( N+J ) = ZERO
                   END IF
                ELSE
                   WORK( J ) = WORK( J )*SQRT( TEMP )
                END IF
             END IF
 30       CONTINUE
  */

  // Some initializations.
  tol3z = sqrt( FLA_Mach_params_opd( FLA_MACH_EPS ) );
  //// printf( "tol3z: %24.17le\n", tol3z );
  //// printf( "eps:   %24.17le\n", FLA_Mach_params_opd( FLA_MACH_EPS ) );
  ptr_d  = buff_d;
  ptr_e  = buff_e;
  ptr_wt = buff_wt;
  ptr_A  = buff_A;

  // Main loop.
  for( j = 0; j < n_A; j++ ) {
    if( * ptr_d != 0.0 ) {
      temp = dabs( * ptr_wt ) / * ptr_d;
      temp = max( 0.0, ( 1.0 + temp ) * ( 1 - temp ) );
      temp5 = * ptr_d / * ptr_e;
      temp2 = temp * temp5 * temp5;
      //// printf( "temp2: %24.17lf   tol3z: %24.17lf\n", temp2, tol3z );
      if( temp2 <= tol3z ) {
        //// printf( "F. Cancel Drm %4d %14.6le %14.6le\n", 
        ////         j, * ptr_wt, * ptr_d );
        if( m_A > 0 ) {
          * ptr_d = dnrm2_( & m_A, ptr_A, & i_one );
          * ptr_e = *ptr_d;
        } else {
          * ptr_d = 0.0;
          * ptr_e = 0.0;
        }
      } else {
        * ptr_d = * ptr_d * sqrt( temp );
      }
    } 
    ptr_A  += ldim_A;
    ptr_d  += st_d;
    ptr_e  += st_e;
    ptr_wt += st_wt;
  }

  return FLA_SUCCESS;
}

/* ========================================================================= */
static int NoFLA_my_idamax( int n, double * buff_data ) {
  int     i, result;
  double  dmax;

  // IDAMAX = 0
  // IF (N.LT.1 .OR. INCX.LE.0) RETURN
  if( n < 1 ) {
    return -1;
  }

  // IDAMAX = 1
  // IF (N.EQ.1) RETURN
  if( n == 1 ) {
    return 0;
  }

  // *
  // * code for increment equal to 1
  // *
  // DMAX = DABS(DX(1))
  // DO I = 2,N
  //   IF (DABS(DX(I)).GT.DMAX) THEN
  //     IDAMAX = I
  //     DMAX = DABS(DX(I))
  //   END IF
  // END DO

  /*
  for( i = 0; i < n; i++ ) {
    printf( "    dabs( buff_data[ %d ] ) = %25.16le\n", 
            i, dabs( buff_data[ i ] ) );
  } */

  result = 0;
  dmax   = dabs( buff_data[ 0 ] );
  for( i = 1; i < n; i++ ) {
    if( dmax < dabs( buff_data[ i ] ) ) {
      result = i;
      dmax   = dabs( buff_data[ i ] );
    }
  }

  return result;
}

































