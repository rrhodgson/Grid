    /*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid 

    Source file: ./lib/algorithms/iterative/MADWF.h

    Copyright (C) 2015

Author: Peter Boyle <paboyle@ph.ed.ac.uk>

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License along
    with this program; if not, write to the Free Software Foundation, Inc.,
    51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.

    See the full license in the file "LICENSE" in the top level distribution directory
    *************************************************************************************/
    /*  END LEGAL */
#pragma once

NAMESPACE_BEGIN(Grid);

template <class Fieldi, class Fieldo,IfNotSame<Fieldi,Fieldo> X=0>
inline void convert(const Fieldi &from,Fieldo &to) 
{
  precisionChange(to,from);
}
template <class Fieldi, class Fieldo,IfSame<Fieldi,Fieldo> X=0>
inline void convert(const Fieldi &from,Fieldo &to) 
{
  to=from;
}

struct MADWFinnerIterCallbackBase{
  virtual void operator()(const RealD current_resid){}
  virtual ~MADWFinnerIterCallbackBase(){}
};

template<class Matrixo,class Matrixi,class PVinverter,class SchurSolver, class Guesser> 
class MADWF 
{
 private:
  typedef typename Matrixo::FermionField FermionFieldo;
  typedef typename Matrixi::FermionField FermionFieldi;

  PVinverter  & PauliVillarsSolvero;// For the outer field
  SchurSolver & SchurSolveri;       // For the inner approx field
  Guesser     & Guesseri;           // To deflate the inner approx solves

  Matrixo & Mato;                   // Action object for outer
  Matrixi & Mati;                   // Action object for inner

  RealD target_resid;
  int   maxiter;

  //operator() is called on "callback" at the end of every inner iteration. This allows for example the adjustment of the inner
  //tolerance to speed up subsequent iteration
  MADWFinnerIterCallbackBase* callback;
  
 public:
  MADWF(Matrixo &_Mato,
	Matrixi &_Mati,
	PVinverter &_PauliVillarsSolvero,
	SchurSolver &_SchurSolveri,
	Guesser & _Guesseri,
	RealD resid,
	int _maxiter,
	MADWFinnerIterCallbackBase* _callback = NULL) :

  Mato(_Mato),Mati(_Mati),
    SchurSolveri(_SchurSolveri),
    PauliVillarsSolvero(_PauliVillarsSolvero),Guesseri(_Guesseri),
    callback(_callback)
    {
      target_resid=resid;
      maxiter     =_maxiter;
    };
   
  void operator() (const FermionFieldo &src4,FermionFieldo &sol5)
  {
    std::cout << GridLogMessage<< " ************************************************" << std::endl;
    std::cout << GridLogMessage<< "  MADWF-like algorithm                           " << std::endl;
    std::cout << GridLogMessage<< " ************************************************" << std::endl;

std::cout << "    FermionFieldi    c0i(Mati.GaugeGrid()); // 4d " << std::endl;
    FermionFieldi    c0i(Mati.GaugeGrid()); // 4d 
std::cout << "    FermionFieldi    y0i(Mati.GaugeGrid()); // 4d" << std::endl;
    FermionFieldi    y0i(Mati.GaugeGrid()); // 4d
std::cout << "    FermionFieldo    c0 (Mato.GaugeGrid()); // 4d " << std::endl;
    FermionFieldo    c0 (Mato.GaugeGrid()); // 4d 
std::cout << "    FermionFieldo    y0 (Mato.GaugeGrid()); // 4d" << std::endl;
    FermionFieldo    y0 (Mato.GaugeGrid()); // 4d

std::cout << "    FermionFieldo    A(Mato.FermionGrid()); // Temporary outer" << std::endl;
    FermionFieldo    A(Mato.FermionGrid()); // Temporary outer
std::cout << "    FermionFieldo    B(Mato.FermionGrid()); // Temporary outer" << std::endl;
    FermionFieldo    B(Mato.FermionGrid()); // Temporary outer
std::cout << "    FermionFieldo    b(Mato.FermionGrid()); // 5d source" << std::endl;
    FermionFieldo    b(Mato.FermionGrid()); // 5d source

std::cout << "    FermionFieldo    c(Mato.FermionGrid()); // PVinv source; reused so store" << std::endl;
    FermionFieldo    c(Mato.FermionGrid()); // PVinv source; reused so store
std::cout << "    FermionFieldo    defect(Mato.FermionGrid()); // 5d source" << std::endl;
    FermionFieldo    defect(Mato.FermionGrid()); // 5d source

std::cout << "    FermionFieldi   ci(Mati.FermionGrid()); " << std::endl;
    FermionFieldi   ci(Mati.FermionGrid()); 
std::cout << "    FermionFieldi   yi(Mati.FermionGrid()); " << std::endl;
    FermionFieldi   yi(Mati.FermionGrid()); 
std::cout << "    FermionFieldi   xi(Mati.FermionGrid()); " << std::endl;
    FermionFieldi   xi(Mati.FermionGrid()); 
std::cout << "    FermionFieldi srci(Mati.FermionGrid()); " << std::endl;
    FermionFieldi srci(Mati.FermionGrid()); 
std::cout << "    FermionFieldi   Ai(Mati.FermionGrid()); " << std::endl;
    FermionFieldi   Ai(Mati.FermionGrid()); 

    RealD m=Mati.Mass();

    ///////////////////////////////////////
    //Import source, include Dminus factors
    ///////////////////////////////////////
std::cout << "    Mato.ImportPhysicalFermionSource(src4,b); " << std::endl;
    Mato.ImportPhysicalFermionSource(src4,b); 
    std::cout << GridLogMessage << " src4 " <<norm2(src4)<<std::endl;
    std::cout << GridLogMessage << " b    " <<norm2(b)<<std::endl;

    defect = b;
    sol5=Zero();
    for (int i=0;i<maxiter;i++) {

      ///////////////////////////////////////
      // Set up c0 from current defect
      ///////////////////////////////////////
std::cout << "      PauliVillarsSolvero(Mato,defect,A);" << std::endl;
      PauliVillarsSolvero(Mato,defect,A);
std::cout << "      Mato.Pdag(A,c);" << std::endl;
      Mato.Pdag(A,c);
std::cout << "      ExtractSlice(c0, c, 0 , 0);" << std::endl;
      ExtractSlice(c0, c, 0 , 0);

      ////////////////////////////////////////////////
      // Solve the inner system with surface term c0
      ////////////////////////////////////////////////
      ci = Zero();  
std::cout << "      convert(c0,c0i); // Possible precison change" << std::endl;
      convert(c0,c0i); // Possible precison change
std::cout << "      InsertSlice(c0i,ci,0, 0);" << std::endl;
      InsertSlice(c0i,ci,0, 0);

      // Dwm P y = Dwm x = D(1) P (c0,0,0,0)^T
std::cout << "      Mati.P(ci,Ai);" << std::endl;
      Mati.P(ci,Ai);
std::cout << "      Mati.SetMass(1.0);      Mati.M(Ai,srci);      Mati.SetMass(m);" << std::endl;
      Mati.SetMass(1.0);      Mati.M(Ai,srci);      Mati.SetMass(m);
std::cout << "      SchurSolveri(Mati,srci,xi,Guesseri); " << std::endl;
      SchurSolveri(Mati,srci,xi,Guesseri); 
std::cout << "      Mati.Pdag(xi,yi);" << std::endl;
      Mati.Pdag(xi,yi);
std::cout << "      ExtractSlice(y0i, yi, 0 , 0);" << std::endl;
      ExtractSlice(y0i, yi, 0 , 0);
std::cout << "      convert(y0i,y0); // Possible precision change" << std::endl;
      convert(y0i,y0); // Possible precision change

      //////////////////////////////////////
      // Propagate solution back to outer system
      // Build Pdag PV^-1 Dm P [-sol4,c2,c3... cL]
      //////////////////////////////////////
      c0 = - y0;
std::cout << "      InsertSlice(c0, c, 0   , 0);" << std::endl;
      InsertSlice(c0, c, 0   , 0);

      /////////////////////////////
      // Reconstruct the bulk solution Pdag PV^-1 Dm P 
      /////////////////////////////
std::cout << "      Mato.P(c,B);" << std::endl;
      Mato.P(c,B);
std::cout << "      Mato.M(B,A);" << std::endl;
      Mato.M(B,A);
std::cout << "      PauliVillarsSolvero(Mato,A,B);" << std::endl;
      PauliVillarsSolvero(Mato,A,B);
std::cout << "      Mato.Pdag(B,A);" << std::endl;
      Mato.Pdag(B,A);

      //////////////////////////////
      // Reinsert surface prop
      //////////////////////////////
std::cout << "      InsertSlice(y0,A,0,0);" << std::endl;
      InsertSlice(y0,A,0,0);

      //////////////////////////////
      // Convert from y back to x 
      //////////////////////////////
std::cout << "      Mato.P(A,B);" << std::endl;
      Mato.P(A,B);

      //         sol5' = sol5 + M^-1 defect
      //               = sol5 + M^-1 src - M^-1 M sol5  ...
      sol5 = sol5 + B;
      std::cout << GridLogMessage << "***************************************" <<std::endl;
      std::cout << GridLogMessage << " Sol5 update "<<std::endl;
      std::cout << GridLogMessage << "***************************************" <<std::endl;
      std::cout << GridLogMessage << " Sol5 now "<<norm2(sol5)<<std::endl;
      std::cout << GridLogMessage << " delta    "<<norm2(B)<<std::endl;

       // New defect  = b - M sol5
std::cout << "       Mato.M(sol5,A);" << std::endl;
       Mato.M(sol5,A);
       defect = b - A;

       std::cout << GridLogMessage << " defect   "<<norm2(defect)<<std::endl;

       double resid = ::sqrt(norm2(defect) / norm2(b));
       std::cout << GridLogMessage << "Residual " << i << ": " << resid  << std::endl;
       std::cout << GridLogMessage << "***************************************" <<std::endl;

       if(callback != NULL) (*callback)(resid);       
       
       if (resid < target_resid) {
	 return;
       }
    }

    std::cout << GridLogMessage << "MADWF : Exceeded maxiter "<<std::endl;
    assert(0);

  }

};

NAMESPACE_END(Grid);
