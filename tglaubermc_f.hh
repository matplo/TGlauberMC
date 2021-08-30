#include "tglaubermc.hh"

namespace TGlauberMC_f
{
//---------------------------------------------------------------------------------
TF1 *getNNProf(Double_t snn=67.6, Double_t omega=0.4, Double_t G=1);

//---------------------------------------------------------------------------------
void runAndSaveNtuple(const Int_t n,
                      const char *sysA        = "Pb",
                      const char *sysB        = "Pb",
                      const Double_t signn    = 67.6,
                      const Double_t sigwidth = -1,
                      const Double_t mind     = 0.4,
		      const Double_t omega    = -1,
                      const Double_t noded    = -1,
                      const char *fname       = 0);

//---------------------------------------------------------------------------------
void runAndSaveNucleons(const Int_t n,                    
                        const char *sysA        = "Pb",           
                        const char *sysB        = "Pb",           
                        const Double_t signn    = 67.6,           
                        const Double_t sigwidth = -1,
                        const Double_t mind     = 0.4,
                        const Bool_t verbose    = 0,
			const Double_t bmin     = 0.0,
			const Double_t bmax     = 20.0,
                        const char *fname       = 0);

//---------------------------------------------------------------------------------
void runAndSmearNtuple(const Int_t n,
                       const Double_t sigs  = 0.4,
                       const char *sysA     = "p",
                       const char *sysB     = "Pb",
                       const Double_t signn = 67.6,
                       const Double_t mind  = 0.4,
		       const Double_t bmin  = 0.0,
		       const Double_t bmax  = 20.0,
                       const char *fname    = 0);


//---------------------------------------------------------------------------------
void runAndOutputLemonTree(const Int_t n,
			   const Double_t sigs  = 0.4,
			   const char *sysA     = "p",
			   const char *sysB     = "Pb",
			   const Double_t signn = 67.6,
			   const Double_t mind  = 0.4,
			   const Double_t bmin  = 0.0,
			   const Double_t bmax  = 20.0,
			   const Bool_t   ogrid = 0,
			   const char *fname    = 0);

//---------------------------------------------------------------------------------
void runAndCalcDens(const Int_t n,
		    const Double_t alpha = 0.1,
		    const char *sysA     = "Pb",
		    const char *sysB     = "Pb",
		    const Double_t signn = 67.6,
		    const Double_t mind  = 0.4,
		    const char *fname    = "glau_dens_hists.root");

};
//---------------------------------------------------------------------------------
