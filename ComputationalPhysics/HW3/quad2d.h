struct NRf2 {
  Doub xsav;
  Doub (*func2d)(const Doub, const Doub);
  Doub operator()(const Doub y)
  {
    return func2d(xsav,y);
  }
};
struct NRf1 {
  Doub (*y1)(Doub);
  Doub (*y2)(Doub);
  NRf2 f2;
  NRf1(Doub yy1(Doub), Doub yy2(Doub)) : y1(yy1),y2(yy2) {}
  Doub operator()(const Doub x)
  {
    f2.xsav=x;
    return qromb(f2,y1(x),y2(x),1.e-9, false);
  }
};

template <class T>
Doub quad2d(T &func, const Doub x1, const Doub x2, Doub y1(Doub), Doub y2(Doub))
{
  NRf1 f1(y1,y2);
  f1.f2.func2d=func;
  Midpnt<NRf1> mid(f1,x1,x2);
  return qromo(mid,1.e-6, false);
}
