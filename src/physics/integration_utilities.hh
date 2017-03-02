#ifndef INTEGRATION_UTILITIES_HH
#define INTEGRATION_UTILITIES_HH

#include <cmath>
#include <vector>
#include <iostream>

//! \brief 1d-adaptive-integration with Gauss 7th Kronrod 15th order

template<typename F>
double integrateG7K15(F &func, double a, double b, double maxerror)
{
  double gauss7 = 0.0000000000000000000000000;
  double kronrod15 = 0.0000000000000000000000000;

  // a has to be the lower bound
  if (a >= b)
    return 0;

  double m = (a+b)/2.;
  double k = (b-a)/2.;

  double funcvalues[15];
  // gauss points are marked with //*
  funcvalues[0] = func.evaluate(0.0000000000000000000000000 + m);       //*
  funcvalues[1] = func.evaluate(0.2077849550078984676006894 * k + m);
  funcvalues[2] = func.evaluate(-0.2077849550078984676006894 * k + m);
  funcvalues[3] = func.evaluate(0.4058451513773971669066064 * k + m);   //*
  funcvalues[4] = func.evaluate(-0.4058451513773971669066064 * k + m);  //*
  funcvalues[5] = func.evaluate(0.5860872354676911302941448 * k + m);
  funcvalues[6] = func.evaluate(-0.5860872354676911302941448 * k + m);
  funcvalues[7] = func.evaluate(0.7415311855993944398638648 * k + m);   //*
  funcvalues[8] = func.evaluate(-0.7415311855993944398638648 * k + m);  //*
  funcvalues[9] = func.evaluate(0.8648644233597690727897128 * k + m);
  funcvalues[10] = func.evaluate(-0.8648644233597690727897128 * k + m);
  funcvalues[11] = func.evaluate(0.9491079123427585245261897 * k + m);  //*
  funcvalues[12] = func.evaluate(-0.9491079123427585245261897 * k + m); //*
  funcvalues[13] = func.evaluate(0.9914553711208126392068547 * k + m);
  funcvalues[14] = func.evaluate(-0.9914553711208126392068547 * k + m);

  // gauss quadrature 7th order
  gauss7 += 0.4179591836734693877551020 * funcvalues[0];
  gauss7 += 0.3818300505051189449503698 * funcvalues[3];
  gauss7 += 0.3818300505051189449503698 * funcvalues[4];
  gauss7 += 0.2797053914892766679014678 * funcvalues[7];
  gauss7 += 0.2797053914892766679014678 * funcvalues[8];
  gauss7 += 0.1294849661688696932706114 * funcvalues[11];
  gauss7 += 0.1294849661688696932706114 * funcvalues[12];
  gauss7 *= k;

  // kronrod quadratur 15th order
  kronrod15 += 0.2094821410847278280129992 * funcvalues[0];
  kronrod15 += 0.2044329400752988924141620 * funcvalues[1];
  kronrod15 += 0.2044329400752988924141620 * funcvalues[2];
  kronrod15 += 0.1903505780647854099132564 * funcvalues[3];
  kronrod15 += 0.1903505780647854099132564 * funcvalues[4];
  kronrod15 += 0.1690047266392679028265834 * funcvalues[5];
  kronrod15 += 0.1690047266392679028265834 * funcvalues[6];
  kronrod15 += 0.1406532597155259187451896 * funcvalues[7];
  kronrod15 += 0.1406532597155259187451896 * funcvalues[8];
  kronrod15 += 0.1047900103222501838398763 * funcvalues[9];
  kronrod15 += 0.1047900103222501838398763 * funcvalues[10];
  kronrod15 += 0.0630920926299785532907007 * funcvalues[11];
  kronrod15 += 0.0630920926299785532907007 * funcvalues[12];
  kronrod15 += 0.0229353220105292249637320 * funcvalues[13];
  kronrod15 += 0.0229353220105292249637320 * funcvalues[14];
  kronrod15 *= k;

  if (std::abs(gauss7 - kronrod15) > maxerror)
    return integrateG7K15(func, a, m, maxerror) + integrateG7K15(func, m, b, maxerror);

  return gauss7;
}

//! \brief 1d-adaptive-integration with Gauss 10th Kronrod 21th order

template<typename F>
double integrateG10K21(F &func,double a, double b, double maxerror)
{
  double gauss10 = 0.0000000000000000000000000;
  double kronrod21 = 0.0000000000000000000000000;

  // a has to be the lower bound
  if (a >= b)
    return 0;

  // scaling from [a,b] to [-1,1]
  double m = (a+b)/2.;
  double k = (b-a)/2.;

  double funcvalues[21];
  // gauss points are marked with //*
  funcvalues[0] = func.evaluate(0.0000000000000000000000000 + m);
  funcvalues[1] = func.evaluate(0.1488743389816312108848260 * k + m);   //*
  funcvalues[2] = func.evaluate(-0.1488743389816312108848260 * k + m);  //*
  funcvalues[3] = func.evaluate(0.2943928627014601981311266 * k + m);
  funcvalues[4] = func.evaluate(-0.2943928627014601981311266 * k + m);
  funcvalues[5] = func.evaluate(0.4333953941292471907992659 * k + m);   //*
  funcvalues[6] = func.evaluate(-0.4333953941292471907992659 * k + m);  //*
  funcvalues[7] = func.evaluate(0.5627571346686046833390001 * k + m);
  funcvalues[8] = func.evaluate(-0.5627571346686046833390001 * k + m);
  funcvalues[9] = func.evaluate(0.6794095682990244062343274 * k + m);   //*
  funcvalues[10] = func.evaluate(-0.6794095682990244062343274 * k + m); //*
  funcvalues[11] = func.evaluate(0.7808177265864168970637176 * k + m);
  funcvalues[12] = func.evaluate(-0.7808177265864168970637176 * k + m);
  funcvalues[13] = func.evaluate(0.8650633666889845107320967 * k + m);  //*
  funcvalues[14] = func.evaluate(-0.8650633666889845107320967 * k + m); //*
  funcvalues[15] = func.evaluate(0.9301574913557082260012072 * k + m);
  funcvalues[16] = func.evaluate(-0.9301574913557082260012072 * k + m);
  funcvalues[17] = func.evaluate(0.9739065285171717200779640 * k + m);  //*
  funcvalues[18] = func.evaluate(-0.9739065285171717200779640 * k + m); //*
  funcvalues[19] = func.evaluate(0.9956571630258080807355273 * k + m);
  funcvalues[20] = func.evaluate(-0.9956571630258080807355273 * k + m);

  // gauss quadrature 10th order
  gauss10 += 0.2955242247147528701738930 * funcvalues[1];
  gauss10 += 0.2955242247147528701738930 * funcvalues[2];
  gauss10 += 0.2692667193099963550912269 * funcvalues[5];
  gauss10 += 0.2692667193099963550912269 * funcvalues[6];
  gauss10 += 0.2190863625159820439955349 * funcvalues[9];
  gauss10 += 0.2190863625159820439955349 * funcvalues[10];
  gauss10 += 0.1494513491505805931457763 * funcvalues[13];
  gauss10 += 0.1494513491505805931457763 * funcvalues[14];
  gauss10 += 0.0666713443086881375935688 * funcvalues[17];
  gauss10 += 0.0666713443086881375935688 * funcvalues[18];
  gauss10 *= k;

  // kronrod quadratur 21th order
  kronrod21 += 0.1494455540029169056649365 * funcvalues[0];
  kronrod21 += 0.1477391049013384913748415 * funcvalues[1];
  kronrod21 += 0.1477391049013384913748415 * funcvalues[2];
  kronrod21 += 0.1427759385770600807970943 * funcvalues[3];
  kronrod21 += 0.1427759385770600807970943 * funcvalues[4];
  kronrod21 += 0.1347092173114733259280540 * funcvalues[5];
  kronrod21 += 0.1347092173114733259280540 * funcvalues[6];
  kronrod21 += 0.1234919762620658510779581 * funcvalues[7];
  kronrod21 += 0.1234919762620658510779581 * funcvalues[8];
  kronrod21 += 0.1093871588022976418992106 * funcvalues[9];
  kronrod21 += 0.1093871588022976418992106 * funcvalues[10];
  kronrod21 += 0.0931254545836976055350655 * funcvalues[11];
  kronrod21 += 0.0931254545836976055350655 * funcvalues[12];
  kronrod21 += 0.0750396748109199527670431 * funcvalues[13];
  kronrod21 += 0.0750396748109199527670431 * funcvalues[14];
  kronrod21 += 0.0547558965743519960313813 * funcvalues[15];
  kronrod21 += 0.0547558965743519960313813 * funcvalues[16];
  kronrod21 += 0.0325581623079647274788190 * funcvalues[17];
  kronrod21 += 0.0325581623079647274788190 * funcvalues[18];
  kronrod21 += 0.0116946388673718742780644 * funcvalues[19];
  kronrod21 += 0.0116946388673718742780644 * funcvalues[20];
  kronrod21 *= k;

  if (std::abs(gauss10 - kronrod21) > maxerror)
    return integrateG10K21(func, a, m, maxerror) + integrateG10K21(func, m, b, maxerror);

  return gauss10;
}

//! \brief 1d-adaptive-integration with Gauss 15th Kronrod 31th order

template<typename F>
double integrateG15K31(F &func,double a, double b, double maxerror)
{
  double gauss15 = 0.0000000000000000000000000;
  double kronrod31 = 0.0000000000000000000000000;

  // a has to be the lower bound
  if (a >= b)
  {
    std::cout << "{integrateG15K31} a >=b" << std::endl;
    return 0;
  }

  double m = (a+b)/2.;
  double k = (b-a)/2.;

  double funcvalues[31];
  // gauss points are marked with //*
  funcvalues[0] = func.evaluate(0.0000000000000000000000000 + m);       //*
  funcvalues[1] = func.evaluate(0.1011420669187174990270742 * k + m);
  funcvalues[2] = func.evaluate(-0.1011420669187174990270742 * k + m);
  funcvalues[3] = func.evaluate(0.2011940939974345223006283 * k + m);   //*
  funcvalues[4] = func.evaluate(-0.2011940939974345223006283 * k + m);  //*
  funcvalues[5] = func.evaluate(0.2991800071531688121667800 * k + m);
  funcvalues[6] = func.evaluate(-0.2991800071531688121667800 * k + m);
  funcvalues[7] = func.evaluate(0.3941513470775633698972074 * k + m);   //*
  funcvalues[8] = func.evaluate(-0.3941513470775633698972074 * k + m);  //*
  funcvalues[9] = func.evaluate(0.4850818636402396806936557 * k + m);
  funcvalues[10] = func.evaluate(-0.4850818636402396806936557 * k + m);
  funcvalues[11] = func.evaluate(0.5709721726085388475372267 * k + m);  //*
  funcvalues[12] = func.evaluate(-0.5709721726085388475372267 * k + m); //*
  funcvalues[13] = func.evaluate(0.6509967412974169705337359 * k + m);
  funcvalues[14] = func.evaluate(-0.6509967412974169705337359 * k + m);
  funcvalues[15] = func.evaluate(0.7244177313601700474161861 * k + m);  //*
  funcvalues[16] = func.evaluate(-0.7244177313601700474161861 * k + m); //*
  funcvalues[17] = func.evaluate(0.7904185014424659329676493 * k + m);
  funcvalues[18] = func.evaluate(-0.7904185014424659329676493 * k + m);
  funcvalues[19] = func.evaluate(0.8482065834104272162006483 * k + m);  //*
  funcvalues[20] = func.evaluate(-0.8482065834104272162006483 * k + m); //*
  funcvalues[21] = func.evaluate(0.8972645323440819008825097 * k + m);
  funcvalues[22] = func.evaluate(-0.8972645323440819008825097 * k + m);
  funcvalues[23] = func.evaluate(0.9372733924007059043077589 * k + m);  //*
  funcvalues[24] = func.evaluate(-0.9372733924007059043077589 * k + m); //*
  funcvalues[25] = func.evaluate(0.9677390756791391342573480 * k + m);
  funcvalues[26] = func.evaluate(-0.9677390756791391342573480 * k + m);
  funcvalues[27] = func.evaluate(0.9879925180204854284895657 * k + m);  //*
  funcvalues[28] = func.evaluate(-0.9879925180204854284895657 * k + m); //*
  funcvalues[29] = func.evaluate(0.9980022986933970602851728 * k + m);
  funcvalues[30] = func.evaluate(-0.9980022986933970602851728 * k + m);

  // gauss quadrature 15th order
  gauss15 += 0.2025782419255612728806202 * funcvalues[0];
  gauss15 += 0.1984314853271115764561183 * funcvalues[1];
  gauss15 += 0.1984314853271115764561183 * funcvalues[2];
  gauss15 += 0.1861610000155622110268006 * funcvalues[3];
  gauss15 += 0.1861610000155622110268006 * funcvalues[4];
  gauss15 += 0.1662692058169939335532009 * funcvalues[5];
  gauss15 += 0.1662692058169939335532009 * funcvalues[6];
  gauss15 += 0.1395706779261543144478048 * funcvalues[7];
  gauss15 += 0.1395706779261543144478048 * funcvalues[8];
  gauss15 += 0.1071592204671719350118695 * funcvalues[9];
  gauss15 += 0.1071592204671719350118695 * funcvalues[10];
  gauss15 += 0.0703660474881081247092674 * funcvalues[11];
  gauss15 += 0.0703660474881081247092674 * funcvalues[12];
  gauss15 += 0.0307532419961172683546284 * funcvalues[13];
  gauss15 += 0.0307532419961172683546284 * funcvalues[14];
  gauss15 *= k;

  // kronrod quadratur 31th order
  kronrod31 += 0.1013300070147915490173748 * funcvalues[0];
  kronrod31 += 0.1007698455238755950449467 * funcvalues[1];
  kronrod31 += 0.1007698455238755950449467 * funcvalues[2];
  kronrod31 += 0.0991735987217919593323932 * funcvalues[3];
  kronrod31 += 0.0991735987217919593323932 * funcvalues[4];
  kronrod31 += 0.0966427269836236785051799 * funcvalues[5];
  kronrod31 += 0.0966427269836236785051799 * funcvalues[6];
  kronrod31 += 0.0931265981708253212254869 * funcvalues[7];
  kronrod31 += 0.0931265981708253212254869 * funcvalues[8];
  kronrod31 += 0.0885644430562117706472754 * funcvalues[9];
  kronrod31 += 0.0885644430562117706472754 * funcvalues[10];
  kronrod31 += 0.0830805028231330210382892 * funcvalues[11];
  kronrod31 += 0.0830805028231330210382892 * funcvalues[12];
  kronrod31 += 0.0768496807577203788944328 * funcvalues[13];
  kronrod31 += 0.0768496807577203788944328 * funcvalues[14];
  kronrod31 += 0.0698541213187282587095201 * funcvalues[15];
  kronrod31 += 0.0698541213187282587095201 * funcvalues[16];
  kronrod31 += 0.0620095678006706402851392 * funcvalues[17];
  kronrod31 += 0.0620095678006706402851392 * funcvalues[18];
  kronrod31 += 0.0534815246909280872653431 * funcvalues[19];
  kronrod31 += 0.0534815246909280872653431 * funcvalues[20];
  kronrod31 += 0.0445897513247648766082273 * funcvalues[21];
  kronrod31 += 0.0445897513247648766082273 * funcvalues[22];
  kronrod31 += 0.0353463607913758462220379 * funcvalues[23];
  kronrod31 += 0.0353463607913758462220379 * funcvalues[24];
  kronrod31 += 0.0254608473267153201868740 * funcvalues[25];
  kronrod31 += 0.0254608473267153201868740 * funcvalues[26];
  kronrod31 += 0.0150079473293161225383748 * funcvalues[27];
  kronrod31 += 0.0150079473293161225383748 * funcvalues[28];
  kronrod31 += 0.0053774798729233489877921 * funcvalues[29];
  kronrod31 += 0.0053774798729233489877921 * funcvalues[30];
  kronrod31 *= k;

  if (std::abs(gauss15 - kronrod31) > maxerror)
    return integrateG15K31(func, a, m, maxerror) + integrateG15K31(func, m, b, maxerror);

  return gauss15;
}

#endif // INTEGRATION_UTILITIES_HH
