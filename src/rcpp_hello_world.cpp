
#include <Rcpp.h>
using namespace Rcpp;

//' The Title! Hello World
//' 
//' The description
//' 
//' @param x Something about a parameter
//' @return Nothing.
//' @details Some details about this function.
// [[Rcpp::export]]
List rcpp_hello_world() {

    CharacterVector x = CharacterVector::create( "foo", "bar" )  ;
    NumericVector y   = NumericVector::create( 0.0, 1.0 ) ;
    List z            = List::create( x, y ) ;
    int blarg = 5;

    return z ;
}
