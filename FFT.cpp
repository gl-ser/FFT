//----------------------------------------------------------------------------//
// Файл FFT.cpp                                                               //
//                                                                            //
//              *** TFFT КЛАСС БЫСТРОГО ПРЕОБРАЗОВАНИЯ ФУРЬЕ ***              //
//                                                                            //
// В терминах обработки сигналов преобразование берёт представление функции   //
// сигнала в виде временных рядов и отображает его в частотный спектр,        //
// где omega — угловая частота.                                               //
// Преобразование превращает функцию времени в функцию частоты.               //
// Преобразование является разложением функции на гармонические составляющие  //
// на различных частотах.                                                     //
// Когда функция f является функцией времени и представляет физический        //
// сигнал, преобразование имеет стандартную интерпретацию как спектр сигнала. //
// Абсолютная величина получающейся в результате комплексной функции F        //
// представляет амплитуды соответствующих частот (omega), в то время как      //
// фазовые сдвиги получаются как аргумент этой комплексной функции.           //
//                                                                            //
// Автор ГЛУЩЕНКО Сергей                                                      //
//                                                                            //
//                                                           Москва, 2022 год //
//----------------------------------------------------------------------------//


#include "FFT.h"


TFFT::TFFT(void)
{
  //Пустой конструктор
}


TFFT::~TFFT(void)
{
  PlusEXP = TEXP();
  MinusEXP = TEXP();
  NN.clear();
}


void TFFT::PreCalcEXP(unsigned int N)
{
  unsigned int temp;
  PlusEXP = TEXP();
  MinusEXP = TEXP();
  NN.clear();

  PlusEXP.resize(static_cast<unsigned int>(floor(log(N)/log(2.0))));
  MinusEXP.resize(static_cast<unsigned int>(floor(log(N)/log(2.0))));
  NN.resize(N+1);

  for(unsigned int i=2; i<NN.size(); i++)
  {
    NN[i] = static_cast<unsigned int>(floor(log(i)/log(2.0)))-1;
  }

  for(unsigned int i=0; i<PlusEXP.size(); i++)
  {
    temp = static_cast<unsigned int>(std::pow(2, i+1));
    PlusEXP[i].n = static_cast<int>(temp);
    MinusEXP[i].n = static_cast<int>(temp);
    PlusEXP[i].data.resize(temp);
    MinusEXP[i].data.resize(temp);

    for(unsigned int wi=0; wi<temp; wi++)
    {
      PlusEXP[i].data[wi] = std::exp(TComplex(0, 2*M_PI*wi/PlusEXP[i].n));
      MinusEXP[i].data[wi] = std::exp(TComplex(0, -2*M_PI*wi/PlusEXP[i].n));
    }
  }
}


TArrComplex TFFT::DirectFFT(const TArrComplex &In)
{
  TArrComplex Out = FFT(In, -1);
  return Out;
}


TArrComplex TFFT::InverseFFT(const TArrComplex &In)
{
  TArrComplex Out = FFT(In, 1);
  double N = static_cast<double>(Out.size());

  for (TArrComplex::iterator itr = Out.begin(); itr != Out.end(); itr++)
  {
    *itr = *itr / N;
  }

  return Out;
}


unsigned int TFFT::MultiplicityOfTwoBig(unsigned int Size)
{
  unsigned int res = static_cast<unsigned int>(floor(log(Size)/log(2.0)));
  return (1<<(res+1));
}


TArrComplex TFFT::FFT(const TArrComplex &In, int Sign)
{
  int i = 0, wi = 0;
  int n = static_cast<int>(In.size());
  TArrComplex A(static_cast<unsigned int>(n)/2), B(static_cast<unsigned int>(n)/2), Out(static_cast<unsigned int>(n));

  if (n == 1)
  {
    return TArrComplex(1, In[0]);
  }

  std::copy_if( In.begin(), In.end(), A.begin(), [&i](TComplex)
  {
    return !(i++ % 2);
  } );
  //В массиве A теперь находится копия четных элементов из массива In

  std::copy_if( In.begin(), In.end(), B.begin(), [&i](TComplex)
  {
    return (i++ % 2);
  } );
  //В массиве B теперь находится копия нечетных элементов из массива In

  TArrComplex At = FFT(A, Sign);
  TArrComplex Bt = FFT(B, Sign);

  for(unsigned int i=0; i<At.size(); i++)
  {
    if (Sign == 1)
      Out[i] = At[i] + Bt[i] * PlusEXP[NN[static_cast<unsigned int>(n)]].data[static_cast<unsigned int>(wi)];
    else
      Out[i] = At[i] + Bt[i] * MinusEXP[NN[static_cast<unsigned int>(n)]].data[static_cast<unsigned int>(wi)];

    wi++;
  }

  for(unsigned int i=0; i<At.size(); i++)
  {
    if (Sign == 1)
      Out[i+static_cast<unsigned int>(n)/2] = At[i] + Bt[i] * PlusEXP[NN[static_cast<unsigned int>(n)]].data[static_cast<unsigned int>(wi)];
    else
      Out[i+static_cast<unsigned int>(n)/2] = At[i] + Bt[i] * MinusEXP[NN[static_cast<unsigned int>(n)]].data[static_cast<unsigned int>(wi)];

    wi++;
  }

  return Out;
}
