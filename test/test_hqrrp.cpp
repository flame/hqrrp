#include <gtest/gtest.h>
#include <hqrrp.h>
#include <util.h>
#include <math.h>

#define RELDTOL 1e-10;
#define ABSDTOL 1e-12;


void expand_pqr(
        double *A,
        int64_t *p,
        double *tau, 
        double *buff_Q,
        int64_t m_A,
        int64_t n_A,
        int64_t ldim_A,
        bool verbose)
{   
    using namespace HQRRP;
    #define max( a, b ) ( (a) > (b) ? (a) : (b) )
    #define min( a, b ) ( (a) < (b) ? (a) : (b) )
    
    int64_t mn_A = min(n_A, m_A);
    // Q is m_A-by-mn_A
    int64_t ldim_Q   = max( 1, m_A );

    lapack::lacpy(
        lapack::MatrixType::General,
        m_A, mn_A, A, ldim_A,
        buff_Q, ldim_Q
    );
    lapack::orgqr(
        m_A, mn_A, mn_A, buff_Q, ldim_Q, tau
    );
    // R is mn_A-by-n_A.
    double *R = (double *) calloc(mn_A * n_A, sizeof(double));
    lapack::lacpy(
        lapack::MatrixType::Upper,
        mn_A, n_A, A, ldim_A,
        R, n_A
    );
    if (verbose) HQRRP::print_double_matrix("R\0", mn_A, n_A, R, mn_A);

    double *temp_A = (double *) calloc(m_A * n_A, sizeof(double));

    blas::gemm(
        blas::Layout::ColMajor,
        blas::Op::NoTrans,
        blas::Op::NoTrans,
        m_A, n_A, mn_A,
        1.0,
        buff_Q, ldim_Q,
        R, mn_A,
        0.0,
        temp_A, m_A
    );

    for (int64_t i = 0; i < n_A; ++i)
    {
        p[i]--;
        int64_t col_idx = p[i];
        blas::copy(m_A, & temp_A[ i* ldim_A] , 1, 
                        & A[ col_idx * ldim_A], 1);
    }
    if (verbose) HQRRP::print_int_vector("p", n_A, p);

    free(temp_A);
    free(R);
}

void test_qrcp(int64_t m_A, int64_t n_A, uint64_t a_seed, bool randalg, bool verbose)
{
    #define max( a, b ) ( (a) > (b) ? (a) : (b) )
    #define min( a, b ) ( (a) < (b) ? (a) : (b) )
    int64_t     nb_alg, pp, info;

    double *buff_A = new double[m_A * n_A];
    double *buff_A_copy = new double[m_A * n_A];
    HQRRP::genmat(m_A, n_A, buff_A, a_seed);
    blas::copy(m_A * n_A, buff_A, 1, buff_A_copy, 1);

    if (verbose) HQRRP::print_double_matrix("mat\0", m_A, n_A, buff_A, m_A);
    
    double mn_A     = min( m_A, n_A );
    int64_t ldim_A   = max( 1, m_A );

    int64_t *buff_p   = ( int64_t * ) calloc( n_A,  sizeof( int64_t ) );
    double *buff_tau = ( double * ) calloc( n_A,  sizeof( double ) );
    
    if (randalg)
    {
        int64_t lwork = max( 1, 128 * n_A );
        double *buff_wk_qp4 = ( double * ) malloc( lwork * sizeof( double ) );
        HQRRP::dgeqp4( & m_A, & n_A, buff_A, & ldim_A, buff_p, buff_tau, 
                        buff_wk_qp4, & lwork, & info );
        if (info != 0) throw blas::Error();
    }
    else
    {
       lapack::geqp3(m_A, n_A, buff_A, ldim_A, buff_p, buff_tau);
    }

    // fill buff_Q, then overwrite buff_A with (Q*R*P').    
    double *buff_Q = (double *) malloc( m_A * mn_A * sizeof(double));
    expand_pqr(buff_A, buff_p, buff_tau, buff_Q, m_A, n_A, ldim_A, verbose);
    if (verbose) HQRRP::print_double_matrix("Q\0", m_A, mn_A, buff_Q, ldim_A);

    // Check that buff_A matches buff_A_copy
    double reldtol = RELDTOL;
    for (int64_t i = 0; i < m_A; ++i)
    {
        for (int64_t j = 0; j < n_A; ++j)
        {
            int64_t ell = i + j*m_A;
            double expect = buff_A_copy[ell];
            double actual = buff_A[ell];
            double abstol = reldtol * min(abs(actual), abs(expect));
            if (abstol == 0.0) abstol = ABSDTOL;
            ASSERT_NEAR(actual, expect, abstol);
        }    
    }

    free(buff_A);
    free(buff_Q);
    free(buff_p);
    free(buff_tau);
    free(buff_A_copy);
}

class TestExpandQRCP : public ::testing::Test 
{
    protected:

    virtual void run(int64_t m_A, int64_t n_A, uint64_t a_seed)
    {
        test_qrcp(m_A, n_A, a_seed, false, false);
    }
};

TEST_F(TestExpandQRCP, Dim10by3)
{
    run(10, 3, 99);
    run(10, 3, 100);
    run(10, 3, 101);
}

TEST_F(TestExpandQRCP, Dim4by4)
{
    run(4, 4, 99);
    run(4, 4, 100);
    run(4, 4, 101);
}

TEST_F(TestExpandQRCP, Dim10by10)
{
    run(10, 10, 99);
    run(10, 10, 100);
    run(10, 10, 101);
}



class TestHQRRP : public ::testing::Test
{
    protected:

    virtual void run(int64_t m_A, int64_t n_A, uint64_t a_seed)
    {
        test_qrcp(m_A, n_A, a_seed, true, false);
    }
};

TEST_F(TestHQRRP, Dim10by3) 
{
   run(10, 3, 99);
   run(10, 3, 100);
   run(10, 3, 101);
}

TEST_F(TestHQRRP, Dim10by10)
{
    run(10, 10, 99);
    run(10, 10, 100);
    run(10, 10, 101);
}

TEST_F(TestHQRRP, Dim300by100)
{
    run(300, 100, 99);
    run(300, 100, 100);
    run(300, 100, 101);
}