#include <gtest/gtest.h>
#include <hqrrp.h>
#include <math.h>

#define RELDTOL 1e-10;
#define ABSDTOL 1e-12;

class TestHQRRP : public ::testing::Test
{
    protected:
        int64_t i;

    virtual void runme(int64_t m, int64_t n)
    {
        i = m;
        // do nothing
    }
};

TEST_F(TestHQRRP, EmptyTest) 
{
    runme(10, 10);
}
