// [[Rcpp::export]]
NumericMatrix abcw(int N, int thin) {

  NumericMatrix res (N,2); // définit une matrice de 0 de taille (N,2) (résultats)
  double x = 0;  // définit un numérique (double précision) initialisé à 0.
  double y = 0;
  int indexline = 0;  // définit un entier qui va nous permettre de remplir la matrice
  for (int i = 0 ; i < N ; i++){
    indexline = i;
    for (int j = 0; j < thin; ++j){
      x = rgamma(1, 3, y * y + 4)[0];
      y = rnorm(1, 1 / (x + 1), 1 / sqrt(2 * ( x + 1)))[0];
    }
    res(indexline,0) = x;
    res(indexline,1) = y;
  }
  return res;
}
