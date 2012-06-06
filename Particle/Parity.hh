#pragma once

class Parity {
public:
  Parity();
  Parity(int parity);
  ~Parity();
  
  int parity();

private:
  int pvalue;
};
