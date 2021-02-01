/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid

Source file: ./tests/solver/Test_zMADWF_prec.cc

Copyright (C) 2015

Author: Christopher Kelly <ckelly@phys.columbia.edu>

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

See the full license in the file "LICENSE" in the top level distribution
directory
*************************************************************************************/
/*  END LEGAL */

//This test computes the zMobius approximation to the Mobius action and uses it within the MADWF context to accelerate an inversion
#include <Grid/Grid.h>

using namespace std;
using namespace Grid;


struct PackRecord
{
    std::string operatorXml, solverXml;
};
struct VecRecord: Serializable
{
    GRID_SERIALIZABLE_CLASS_MEMBERS(VecRecord,
                                    unsigned int, index,
                                    double,       eval);
    VecRecord(void): index(0), eval(0.) {}
};

inline void readHeader(PackRecord &record, ScidacReader &binReader)
    {
        std::string recordXml;

        binReader.readLimeObject(recordXml, SCIDAC_FILE_XML);
        XmlReader xmlReader(recordXml, true, "eigenPackPar");
        xmlReader.push();
        xmlReader.readCurrentSubtree(record.operatorXml);
        xmlReader.nextElement();
        xmlReader.readCurrentSubtree(record.solverXml);
    }
template <typename T, typename TIo = T>
void readElement(T &evec, RealD &eval, const unsigned int index,
                 ScidacReader &binReader, TIo *ioBuf = nullptr)
{
    VecRecord vecRecord;

    std::cout << "Reading eigenvector " << index << std::endl;
    if (ioBuf == nullptr)
    {
      std::cout << "ioBuf == nullptr" << std::endl;
        binReader.readScidacFieldRecord(evec, vecRecord);
      std::cout << "readScidacFieldRecord complete" << std::endl;
    }
    else
    {
      std::cout << "ioBuf != nullptr" << std::endl;
        binReader.readScidacFieldRecord(*ioBuf, vecRecord);
      std::cout << "readScidacFieldRecord complete" << std::endl;
        precisionChange(evec, *ioBuf);
      std::cout << "prec change complete" << std::endl;
    }
    if (vecRecord.index != index)
    {
        std::cout << "Eigenvector " + std::to_string(index) + " has a"
                        + " wrong index (expected " + std::to_string(vecRecord.index) 
                        + ")" << std::endl;
    }
    eval = vecRecord.eval;
      std::cout << "got eval" << std::endl;
}

template <typename T, typename TIo = T>
static void readPack(std::vector<T> &evec, std::vector<RealD> &eval,
                     PackRecord &record, const std::string filename, 
                     const unsigned int size, bool multiFile, 
                     GridBase *gridIo = nullptr)
{
    std::unique_ptr<TIo> ioBuf{nullptr};
    ScidacReader         binReader;

    if (std::is_same<T,TIo>::value == false)
    {
        std::cout << "load type != use type" << std::endl;
        if (gridIo == nullptr)
        {
            std::cout << "I/O type different from vector type but null I/O grid passed" << std::endl;
        }
        ioBuf.reset(new TIo(gridIo));
    }

    if (multiFile)
    {
      std::cout << "multiFile" << std::endl;
        std::string fullFilename;

        for(int k = 0; k < size; ++k) 
        {
            fullFilename = filename + "/v" + std::to_string(k) + ".bin";
            binReader.open(fullFilename);
            readHeader(record, binReader);
            readElement(evec[k], eval[k], k, binReader, ioBuf.get());
            binReader.close();
        }
    }
    else
    {
        binReader.open(filename);
        readHeader(record, binReader);
        for(int k = 0; k < size; ++k) 
        {
            readElement(evec[k], eval[k], k, binReader, ioBuf.get());
        }
        binReader.close();
    }
}




struct TestParams{
  bool load_config;
  std::string config_file;

  std::string evec_file;
  int esize;
  bool multiFile;


  double mass;

  std::string outer_precon;
  std::string inner_precon;
  
  int Ls_outer;
  double b_plus_c_outer;
  double resid_outer;
  int itter_outer;
  
  int Ls_inner;
  double b_plus_c_inner; //irrelevant for ZMobius
  double resid_inner;
  int itter_inner;
  bool zmobius_inner;
  bool accept_gammas;
  double lambda_max; //upper bound of H_T eigenvalue range required to generate zMobius approximation

  double boundary_t;

  
  TestParams(): load_config(true), config_file("ckpoint_lat.1000"), mass(0.01),
    Ls_outer(24), b_plus_c_outer(2.0), resid_outer(1e-8), itter_outer(100), 
    Ls_inner(10), b_plus_c_inner(1.0), resid_inner(1e-8), itter_inner(30000), 
    zmobius_inner(true), accept_gammas(false), lambda_max(1.42), 
    outer_precon("Standard"), inner_precon("Standard"), evec_file(""), esize(1), multiFile(false), boundary_t(-1.)
  {}
  
  void write(const std::string &file) const{
    XmlWriter wr(file);
#define DOIT(A) wr.writeDefault(#A, A)
    DOIT(load_config);
    DOIT(config_file);
    DOIT(mass);
    DOIT(outer_precon);
    DOIT(inner_precon);
    DOIT(Ls_outer);
    DOIT(b_plus_c_outer);
    DOIT(resid_outer);
    DOIT(itter_outer);
    DOIT(Ls_inner);
    DOIT(b_plus_c_inner);
    DOIT(resid_inner);
    DOIT(itter_inner);
    DOIT(zmobius_inner);
    DOIT(accept_gammas);
    DOIT(lambda_max);
    DOIT(evec_file);
    DOIT(esize);
    DOIT(multiFile);
    DOIT(boundary_t);
#undef DOIT
  }
  void read(const std::string &file){
    XmlReader rd(file);
#define DOIT(A) rd.readDefault(#A, A)
    DOIT(load_config);
    DOIT(config_file);
    DOIT(mass);
    DOIT(outer_precon);
    DOIT(inner_precon);
    DOIT(Ls_outer);
    DOIT(b_plus_c_outer);
    DOIT(resid_outer);
    DOIT(Ls_inner);
    DOIT(b_plus_c_inner);
    DOIT(resid_inner);
    DOIT(zmobius_inner);
    DOIT(accept_gammas);
    DOIT(lambda_max);
    DOIT(evec_file);
    DOIT(esize);
    DOIT(multiFile);
    DOIT(boundary_t);
#undef DOIT
  }
};

struct RunParamsPrecStd{
  typedef SchurRedBlackDiagMooeeSolve<LatticeFermionD> SchurSolverType;

  template<typename Action>
  using HermOpType = SchurDiagMooeeOperator<Action, LatticeFermionD>;
};

struct RunParamsPrecDiagTwo{
  typedef SchurRedBlackDiagTwoSolve<LatticeFermionD> SchurSolverType;

  template<typename Action>
  using HermOpType = SchurDiagTwoOperator<Action, LatticeFermionD>;
};


struct CGincreaseTol : public MADWFinnerIterCallbackBase{
  ConjugateGradient<LatticeFermionD> &cg_inner;  
  RealD outer_resid;

  CGincreaseTol(ConjugateGradient<LatticeFermionD> &cg_inner,
    RealD outer_resid): cg_inner(cg_inner), outer_resid(outer_resid){}
  
  void operator()(const RealD current_resid){
    std::cout << "CGincreaseTol with current residual " << current_resid << " changing inner tolerance " << cg_inner.Tolerance << " -> ";
    while(cg_inner.Tolerance < current_resid) cg_inner.Tolerance *= 2;    
    //cg_inner.Tolerance = outer_resid/current_resid;
    std::cout << cg_inner.Tolerance << std::endl;
  }
};









template<typename RunParamsOuter, typename RunParamsInner>
void run(const TestParams &params){

  std::cout << "Loading config file '" << params.config_file << "'" << std::endl
            << "Loading evec file '" << params.evec_file << "'" << std::endl
            << "Using    outer Ls = " << params.Ls_outer << std::endl
            << "Using    inner Ls = " << params.Ls_inner << std::endl
            << "Using        mass = " << params.mass << std::endl
            << "Using itter inner = " << params.itter_inner << std::endl
            << "Using itter outer = " << params.itter_outer << std::endl;


  RealD bmc = 1.0; //use Shamir kernel
  std::vector<ComplexD> gamma_inner;

  if (params.accept_gammas) {
    std::cout << "Accept parameters" << std::endl;
    // gamma_inner.push_back(ComplexD(1.458064389850479e+00,-0.000000000000000e+00));
    // gamma_inner.push_back(ComplexD(1.182313183893475e+00,-0.000000000000000e+00));
    // gamma_inner.push_back(ComplexD(8.309511666859551e-01,-0.000000000000000e+00));
    // gamma_inner.push_back(ComplexD(5.423524091567911e-01,-0.000000000000000e+00));
    // gamma_inner.push_back(ComplexD(3.419850204537295e-01,-0.000000000000000e+00));
    // gamma_inner.push_back(ComplexD(2.113790261902896e-01,-0.000000000000000e+00));
    // gamma_inner.push_back(ComplexD(1.260742995029118e-01,-0.000000000000000e+00));
    // gamma_inner.push_back(ComplexD(9.901366519626265e-02,-0.000000000000000e+00));
    // gamma_inner.push_back(ComplexD(6.863249884465925e-02, 5.506585308274019e-02));
    // gamma_inner.push_back(ComplexD(6.863249884465925e-02,-5.506585308274019e-02));


    gamma_inner.push_back(ComplexD(1.45806438985048,0));
    gamma_inner.push_back(ComplexD(1.18231318389348,0));
    gamma_inner.push_back(ComplexD(0.830951166685955,0));
    gamma_inner.push_back(ComplexD(0.542352409156791,0));
    gamma_inner.push_back(ComplexD(0.341985020453729,0));
    gamma_inner.push_back(ComplexD(0.21137902619029,0));
    gamma_inner.push_back(ComplexD(0.126074299502912,0));
    gamma_inner.push_back(ComplexD(0.0990136651962626,0)); 
    gamma_inner.push_back(ComplexD(0.0686324988446592,0.0550658530827402));
    gamma_inner.push_back(ComplexD(0.0686324988446592,-0.0550658530827402));

  } else {
    std::cout << "Compute parameters" << std::endl;
    if(params.zmobius_inner){
      Approx::computeZmobiusGamma(gamma_inner, params.b_plus_c_inner, params.Ls_inner, params.b_plus_c_outer, params.Ls_outer, params.lambda_max);
    }else{
      Approx::zolotarev_data *zdata = Approx::higham(1.0,params.Ls_inner);
      gamma_inner.resize(params.Ls_inner);
      for(int s=0;s<params.Ls_inner;s++) gamma_inner[s] = zdata->gamma[s];
      Approx::zolotarev_free(zdata);
    }
  }

  std::cout << "gamma:\n";
  for(int s=0;s<params.Ls_inner;s++) std::cout << s << " " << gamma_inner[s] << std::endl;


  std::cout << GridLogMessage << "Lattice dimensions: " << GridDefaultLatt()
            << "   Ls: " << params.Ls_outer << std::endl;

  RealD M5 = 1.8;

  RealD b_outer = (params.b_plus_c_outer + bmc)/2.;
  RealD c_outer = (params.b_plus_c_outer - bmc)/2.;

  RealD b_inner = (params.b_plus_c_inner + bmc)/2.;
  RealD c_inner = (params.b_plus_c_inner - bmc)/2.;






  // Create Grids
  GridCartesian* UGrid = SpaceTimeGrid::makeFourDimGrid(
      GridDefaultLatt(), GridDefaultSimd(Nd, vComplexD::Nsimd()),
      GridDefaultMpi());
  GridRedBlackCartesian* UrbGrid = SpaceTimeGrid::makeFourDimRedBlackGrid(UGrid);


  GridCartesian* FGrid_outer = SpaceTimeGrid::makeFiveDimGrid(params.Ls_outer, UGrid);
  GridCartesian* FGrid_inner = SpaceTimeGrid::makeFiveDimGrid(params.Ls_inner, UGrid);

  GridRedBlackCartesian* FrbGrid_outer = SpaceTimeGrid::makeFiveDimRedBlackGrid(params.Ls_outer, UGrid);
  GridRedBlackCartesian* FrbGrid_inner = SpaceTimeGrid::makeFiveDimRedBlackGrid(params.Ls_inner, UGrid);


  GridRedBlackCartesian* gridIo = SpaceTimeGrid::makeFiveDimRedBlackGrid(params.Ls_inner, UGrid);


  // Create source
  std::vector<int> seeds4({1, 2, 3, 4});

  GridParallelRNG RNG4(UGrid);
  RNG4.SeedFixedIntegers(seeds4);

  LatticeFermionD src4(UGrid); 
  src4 = Zero();
  random(RNG4,src4);
  // Coordinate c(std::vector<int>({0,0,0,0}));
  // iScalar<iVector<iVector<ComplexD, 3>, 4> > s;
  // s = 1.;
  // pokeSite(s, src4, c);



  // Get/make Gauge field
  LatticeGaugeFieldD Umu(UGrid);
  if(params.load_config){
    FieldMetaData header;
    NerscIO::readConfiguration(Umu, header, params.config_file);

    for(int i=0;i<Nd;i++){
      assert(header.dimension[i] == GridDefaultLatt()[i]);
    }
  }else{    
    SU<Nc>::HotConfiguration(RNG4, Umu);
  }
    










  WilsonImplParams implParams;
  implParams.boundary_phases[3] = params.boundary_t;
  std::cout << implParams.boundary_phases << std::endl;

  MobiusFermionD  D_outer(Umu, *FGrid_outer, *FrbGrid_outer, *UGrid, *UrbGrid, params.mass, M5,              b_outer, c_outer, implParams);
  ZMobiusFermionD D_inner(Umu, *FGrid_inner, *FrbGrid_inner, *UGrid, *UrbGrid, params.mass, M5, gamma_inner, b_inner, c_inner, implParams);

  LatticeFermionD src_outer(FGrid_outer);
  D_outer.ImportPhysicalFermionSource(src4,src_outer); //applies D_- 

  LatticeFermionD result_Mob(FGrid_outer);
  LatticeFermionD result_Mob_MADWF(FGrid_outer);
  LatticeFermionD result_diff(FGrid_outer);


  GridStopWatch CGTimer;
  GridStopWatch CGTimer_MADWF;

{
  // Regular Mobuis solve

  // Setup outer correction CG
  ConjugateGradient<LatticeFermionD> CG_outer_correction(params.resid_outer, params.itter_inner);
  typename RunParamsOuter::SchurSolverType SchurSolver_outer(CG_outer_correction, false, true);

  ZeroGuesser<LatticeFermion> zeroGuess;

  result_Mob = Zero();

    CGTimer.Start();
  SchurSolver_outer(D_outer, src_outer, result_Mob, zeroGuess);

    CGTimer.Stop();
}

{
  // MADWF

  // Create Outer PV_CG
  ConjugateGradient<LatticeFermionD> CG_outer(params.resid_outer, params.itter_inner);
  // typedef PauliVillarsSolverFourierAccel<LatticeFermionD, LatticeGaugeFieldD> PVtype;
  // PVtype PV_outer(Umu, CG_outer);

  // typedef PauliVillarsSolverUnprec<LatticeFermionD> PVtype;
  // PVtype PV_outer(CG_outer);
  typename RunParamsOuter::SchurSolverType PV_SchurSolver_outer(CG_outer);
  typedef PauliVillarsSolverRBprec<LatticeFermionD,typename RunParamsOuter::SchurSolverType> PVtype;
  PVtype PV_outer(PV_SchurSolver_outer);

  // Create Inner Shur_CG
  ConjugateGradient<LatticeFermionD> CG_inner(params.resid_inner, params.itter_inner, 0);
  typename RunParamsInner::SchurSolverType SchurSolver_inner(CG_inner);

  CGincreaseTol update(CG_inner, params.resid_outer);


  // Create Results container
  result_Mob_MADWF = Zero();


  if (params.evec_file.empty()) {
    // Zero Guesser 
    ZeroGuesser<LatticeFermion> zeroGuess;

    // Create MADWF
    MADWF<MobiusFermionD, ZMobiusFermionD, PVtype, typename RunParamsInner::SchurSolverType, ZeroGuesser<LatticeFermion> > 
        madwf(D_outer, D_inner, PV_outer, SchurSolver_inner, zeroGuess, params.resid_outer, params.itter_outer, &update);


    // Run MADWF
    CGTimer_MADWF.Start();
    madwf(src_outer, result_Mob_MADWF);
    // madwf(src4, result_Mob_MADWF);
    CGTimer_MADWF.Stop();

  } else {
    // Deflated Guesser
    std::vector<WilsonImplD::FermionField> evec(params.esize, FrbGrid_inner);
    std::vector<RealD> eval(params.esize);
    PackRecord         record;

    readPack<WilsonImplD::FermionField, WilsonImplD::FermionField>(evec, eval,
                       record, params.evec_file, 
                       params.esize, params.multiFile,
                       gridIo);
    DeflatedGuesser<LatticeFermion> difGuess(evec, eval);


    // Create MADWF
    MADWF<MobiusFermionD, ZMobiusFermionD, PVtype, typename RunParamsInner::SchurSolverType, DeflatedGuesser<LatticeFermion> > 
        madwf(D_outer, D_inner, PV_outer, SchurSolver_inner, difGuess,  params.resid_outer, params.itter_outer, &update);


    // Run MADWF
    CGTimer_MADWF.Start();
    madwf(src_outer, result_Mob_MADWF);
    // madwf(src4, result_Mob_MADWF);
    CGTimer_MADWF.Stop();
  }


  result_diff = result_Mob - result_Mob_MADWF;
  std::cout << GridLogMessage<< "MADWF (pre correction)  result norm " << norm2(result_Mob_MADWF) << std::endl
            << GridLogMessage<< "Mobius result norm " << norm2(result_Mob) << std::endl
            << GridLogMessage<< "diff result norm   " << norm2(result_diff) << std::endl
            << std::endl;


  // Setup outer correction CG
  ConjugateGradient<LatticeFermionD> CG_outer_correction(params.resid_outer, params.itter_inner);
  typename RunParamsOuter::SchurSolverType SchurSolver_outer(CG_outer, false, true);


  SchurSolver_outer(D_outer, src_outer, result_Mob_MADWF);
}

  result_diff = result_Mob - result_Mob_MADWF;


  std::cout << GridLogMessage << "Total MADWF time  : " << CGTimer_MADWF.Elapsed()
            << std::endl;
  std::cout << GridLogMessage << "Total Mobius time : " << CGTimer.Elapsed()
            << std::endl;

  std::cout << GridLogMessage<< "MADWF  result norm " << norm2(result_Mob_MADWF) << std::endl
            << GridLogMessage<< "Mobius result norm " << norm2(result_Mob) << std::endl
            << GridLogMessage<< "diff result norm   " << norm2(result_diff) << std::endl
      << "" << std::endl;


  // std::cout << GridLogMessage << "######## Dhop calls summary : outer" << std::endl;
  // D_outer.Report();
  // std::cout << GridLogMessage << "######## Dhop calls summary : inner" << std::endl;
  // D_inner.Report();


}



int main(int argc, char** argv) {
  std::cout << "Init" << std::endl;
  Grid_init(&argc, &argv);

  TestParams params;
  
  if( GridCmdOptionExists(argv,argv+argc,"--params") ){
    std::string pfile = GridCmdOptionPayload(argv,argv+argc,"--params");
    if(pfile == "TEMPLATE"){
      params.write("params.templ");
      return 0;
    }else{
      params.read(pfile);
    }
  }

  // if(params.outer_precon == "Standard" && params.inner_precon == "Standard" ){
  //   run<RunParamsPrecStd, RunParamsPrecStd>(params);
  // }else if(params.outer_precon == "DiagTwo" && params.inner_precon == "Standard"){
  //   run<RunParamsPrecDiagTwo, RunParamsPrecStd>(params);
  // }else if(params.outer_precon == "Standard" && params.inner_precon == "DiagTwo"){
  //   run<RunParamsPrecStd, RunParamsPrecDiagTwo>(params);
  // }else if(params.outer_precon == "DiagTwo" && params.inner_precon == "DiagTwo"){
    run<RunParamsPrecDiagTwo, RunParamsPrecDiagTwo>(params);
  // }else assert(0);

  Grid_finalize();
}
