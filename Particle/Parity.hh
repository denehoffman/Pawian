#ifndef PARITY_HH
#define PARITY_HH

class Parity {

public:
  Parity();
  Parity(int parity);
  ~Parity();

  int parity();

private:
  int pvalue;

};

#endif
