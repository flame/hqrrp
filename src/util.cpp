#include <util.h>
#include <Random123/philox.h>
#include <Random123/uniform.hpp>

namespace HQRRP {

void genmat(
	int64_t n_rows,
	int64_t n_cols,
	double* mat,
	uint64_t seed)
{
	typedef r123::Philox2x64 CBRNG;
	CBRNG::key_type key = {{seed}};
	CBRNG::ctr_type ctr = {{0,0}};
	CBRNG g;
	uint64_t prod = n_rows * n_cols;
	for (uint64_t i = 0; i < prod; ++i)
	{
		ctr[0] = i;
		CBRNG::ctr_type rand = g(ctr, key);
		mat[i] = r123::uneg11<double>(rand.v[0]);
	}
}

// ============================================================================
void print_double_matrix( char * name, int64_t m_A, int64_t n_A, 
                double * buff_A, int64_t ldim_A )
{
  int64_t  i, j;

  printf( "%s = np.array([\n", name );
  for( i = 0; i < m_A; i++ )
  {
    printf(" [");
    for( j = 0; j < n_A - 1; j++ )
	  {
      printf( "%le, ", buff_A[ i + j * ldim_A ] );
    }
    // ++j;
    printf("%le ],\n", buff_A[ i + j * ldim_A ]);
  }
  printf( "])\n" );
}

// ============================================================================
void print_double_vector( char * name, int64_t n_v, double * buff_v )
{
  int64_t  i, j;

  printf( "%s = [\n", name );
  for( i = 0; i < n_v; i++ )
  {
    printf( "%le\n", buff_v[ i ] );
  }
  printf( "\n" );
  printf( "];\n" );
}

// ============================================================================
void print_int_vector( char * name, int64_t n_v, int64_t * buff_v )
{
  int64_t  i, j;

  printf( "%s = [\n", name );
  for( i = 0; i < n_v; i++ )
  {
    printf( "%d\n", (int) buff_v[ i ] );
  }
  printf( "];\n" );
}

// ============================================================================
void init_pvt( int64_t n_p, int64_t * buff_p )
{
  int64_t  i;

  for( i = 0; i < n_p; i++ )
  {
    buff_p[ i ] = ( i + 1 );
  }
}

// ============================================================================
void set_pvt_to_zero( int64_t n_p, int64_t * buff_p )
{
  int64_t  i;

  for( i = 0; i < n_p; i++ )
  {
    buff_p[ i ] = 0;
  }
}


}  // end namespace HQRRCP
