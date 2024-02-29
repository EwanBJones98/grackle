#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <unistd.h>
#include <hdf5.h>

#include <read_simba_snapshot.h>

// This is copied directly from gizmo. It is needed to read the snapshots.
void get_dataset_name(enum iofields blocknr, char *buf)
{
    strcpy(buf, "default");
    
    switch (blocknr)
    {
        case IO_POS:
            strcpy(buf, "Coordinates");
            break;
        case IO_VEL:
            strcpy(buf, "Velocities");
            break;
        case IO_ID:
            strcpy(buf, "ParticleIDs");
            break;
        case IO_MASS:
            strcpy(buf, "Masses");
            break;
        case IO_DUSTMASS:
            strcpy(buf, "Dust_Masses");
            break;
        case IO_U:
            strcpy(buf, "InternalEnergy");
            break;
        case IO_UISM:
            strcpy(buf, "InternalEnergyISM");
            break;
        case IO_RHO:
            strcpy(buf, "Density");
            break;
        case IO_RHOISM:
            strcpy(buf, "DensityISM");
            break;
        case IO_NE:
            strcpy(buf, "ElectronAbundance");
            break;
        case IO_NH:
            strcpy(buf, "NeutralHydrogenAbundance");
            break;
        case IO_TMAX:
            strcpy(buf, "TemperatureMax");
            break;
        case IO_RADGAMMA:
            strcpy(buf, "PhotonEnergy");
            break;
        case IO_RAD_ACCEL:
            strcpy(buf, "RadiativeAcceleration");
            break;
        case IO_HII:
            strcpy(buf, "HII");
            break;
        case IO_HeI:
            strcpy(buf, "HeI");
            break;
        case IO_HeII:
            strcpy(buf, "HeII");
            break;
        case IO_HeIII:
            strcpy(buf, "HeIII");
            break;
        case IO_H2I:
            strcpy(buf, "H2I");
            break;
        case IO_H2II:
            strcpy(buf, "H2II");
            break;
        case IO_HM:
            strcpy(buf, "HM");
            break;
        case IO_HD:
            strcpy(buf, "HD  ");
            break;
        case IO_DI:
            strcpy(buf, "DI  ");
            break;
        case IO_DII:
            strcpy(buf, "DII ");
            break;
        case IO_HeHII:
            strcpy(buf, "HeHp");
            break;
        case IO_DELAYTIME:
            strcpy(buf, "DelayTime");
            break;
        case IO_HSML:
            strcpy(buf, "SmoothingLength");
            break;
        case IO_SFR:
            strcpy(buf, "StarFormationRate");
            break;
        case IO_AGE:
            strcpy(buf, "StellarFormationTime");
            break;
        case IO_INITMASS:
            strcpy(buf, "StellarInitMass");
            break;
        case IO_GRAINSIZE:
            strcpy(buf, "GrainSize");
            break;
        case IO_HSMS:
            strcpy(buf, "StellarSmoothingLength");
            break;
        case IO_Z:
            strcpy(buf, "Metallicity");
            break;
		case IO_DUSTZ:
			strcpy(buf, "Dust_Metallicity");
			break;
        case IO_POT:
            strcpy(buf, "Potential");
            break;
        case IO_ACCEL:
            strcpy(buf, "Acceleration");
            break;
        case IO_DTENTR:
            strcpy(buf, "RateOfChangeOfInternalEnergy");
            break;
        case IO_STRESSDIAG:
            strcpy(buf, "DiagonalStressTensor");
            break;
        case IO_STRESSOFFDIAG:
            strcpy(buf, "OffDiagonalStressTensor");
            break;
        case IO_STRESSBULK:
            strcpy(buf, "BulkStressTensor");
            break;
        case IO_SHEARCOEFF:
            strcpy(buf, "ShearCoefficient");
            break;
        case IO_TSTP:
            strcpy(buf, "TimeStep");
            break;
        case IO_BFLD:
            strcpy(buf, "MagneticField");
            break;
        case IO_DBDT:
            strcpy(buf, "RateOfChangeOfMagneticField");
            break;
        case IO_VRMS:
            strcpy(buf, "RMSVelocity");
            break;
        case IO_VBULK:
            strcpy(buf, "BulkVelocity");
            break;
        case IO_VRAD:
            strcpy(buf, "RMSRadialVelocity");
            break;
        case IO_VTAN:
            strcpy(buf, "RMSTangentialVelocity");
            break;
        case IO_TRUENGB:
            strcpy(buf, "TrueNumberOfNeighbours");
            break;
        case IO_VDIV:
            strcpy(buf, "VelocityDivergence");
            break;
        case IO_VROT:
            strcpy(buf, "VelocityCurl");
            break;
        case IO_VORT:
            strcpy(buf, "Vorticity");
            break;
        case IO_IMF:
            strcpy(buf, "IMFFormationProperties");
            break;
        case IO_COSMICRAY_ENERGY:
            strcpy(buf, "CosmicRayEnergy");
            break;
        case IO_DIVB:
            strcpy(buf, "DivergenceOfMagneticField");
            break;
        case IO_ABVC:
            strcpy(buf, "ArtificialViscosity");
            break;
        case IO_AMDC:
            strcpy(buf, "ArtMagneticDissipation");
            break;
        case IO_PHI:
            strcpy(buf, "DivBcleaningFunctionPhi");
            break;
        case IO_GRADPHI:
            strcpy(buf, "DivBcleaningFunctionGradPhi");
            break;
        case IO_ROTB:
            strcpy(buf, "RotationB");
            break;
        case IO_COOLRATE:
            strcpy(buf, "CoolingRate");
            break;
        case IO_CONDRATE:
            strcpy(buf, "ConductionRate");
            break;
        case IO_DENN:
            strcpy(buf, "Denn");
            break;
        case IO_BHMASS:
            strcpy(buf, "BH_Mass");
            break;
        case IO_BH_DIST:
            strcpy(buf, "BH_Dist");
            break;
        case IO_BHMASSALPHA:
            strcpy(buf, "BH_Mass_AlphaDisk");
            break;
        case IO_ACRB:
            strcpy(buf, "BH_AccretionLength");
            break;
        case IO_BHMDOT:
            strcpy(buf, "BH_Mdot");
            break;
        case IO_BHPROGS:
            strcpy(buf, "BH_NProgs");
            break;
        case IO_BHMBUB:
            strcpy(buf, "BH_Mass_bubbles");
            break;
        case IO_BHMINI:
            strcpy(buf, "BH_Mass_ini");
            break;
        case IO_BHMRAD:
            strcpy(buf, "BH_Mass_radio");
            break;
        case IO_TIDALTENSORPS:
            strcpy(buf, "TidalTensorPS");
            break;
        case IO_DISTORTIONTENSORPS:
            strcpy(buf, "DistortionTensorPS");
            break;
        case IO_CAUSTIC_COUNTER:
            strcpy(buf, "CausticCounter");
            break;
        case IO_FLOW_DETERMINANT:
            strcpy(buf, "FlowDeterminant");
            break;
        case IO_STREAM_DENSITY:
            strcpy(buf, "StreamDensity");
            break;
        case IO_SECONDORDERMASS:
            strcpy(buf, "2lpt-mass");
            break;
        case IO_PHASE_SPACE_DETERMINANT:
            strcpy(buf, "PhaseSpaceDensity");
            break;
        case IO_ANNIHILATION_RADIATION:
            strcpy(buf, "AnnihilationRadiation");
            break;
        case IO_LAST_CAUSTIC:
            strcpy(buf, "LastCaustic");
            break;
        case IO_SHEET_ORIENTATION:
            strcpy(buf, "SheetOrientation");
            break;
        case IO_INIT_DENSITY:
            strcpy(buf, "InitDensity");
            break;
        case IO_EOSTEMP:
            strcpy(buf, "Temperature");
            break;
        case IO_EOSABAR:
            strcpy(buf, "Abar");
            break;
        case IO_EOSYE:
            strcpy(buf, "Ye");
            break;
        case IO_PRESSURE:
            strcpy(buf, "Pressure");
            break;
        case IO_EDDINGTON_TENSOR:
            strcpy(buf, "EddingtonTensor");
            break;
        case IO_DMHSML:
            strcpy(buf, "DM Hsml");
            break;
        case IO_DMDENSITY:
            strcpy(buf, "DM Density");
            break;
        case IO_DMVELDISP:
            strcpy(buf, "DM Velocity Dispersion");
            break;
        case IO_DMHSML_V:
            strcpy(buf, "DM Hsml Voronoi");
            break;
        case IO_DMDENSITY_V:
            strcpy(buf, "DM Density Voronoi");
            break;
        case IO_CHEM:
            strcpy(buf, "ChemicalAbundances");
            break;
        case IO_AGS_SOFT:
            strcpy(buf, "AGS-Softening");
            break;
        case IO_AGS_ZETA:
            strcpy(buf, "AGS-Zeta");
            break;
        case IO_AGS_OMEGA:
            strcpy(buf, "AGS-Omega");
            break;
        case IO_AGS_CORR:
            strcpy(buf, "AGS-Correction");
            break;
        case IO_AGS_NGBS:
            strcpy(buf, "AGS-Neighbours");
            break;
        case IO_VSTURB_DISS:
            strcpy(buf, "TurbulenceDissipation");
            break;
        case IO_VSTURB_DRIVE:
            strcpy(buf, "TurbulenceDriving");
            break;
        case IO_MG_PHI:
            strcpy(buf, "ModifiedGravityPhi");
            break;
        case IO_MG_ACCEL:
            strcpy(buf, "ModifiedGravityAcceleration");
            break;
        case IO_grHI:
            strcpy(buf, "GrackleHI");
            break;
        case IO_grHII:
            strcpy(buf, "GrackleHII");
            break;
        case IO_grHM:
            strcpy(buf, "GrackleHM");
            break;
        case IO_grHeI:
            strcpy(buf, "GrackleHeI");
            break;
        case IO_grHeII:
            strcpy(buf, "GrackleHeII");
            break;
        case IO_grHeIII:
            strcpy(buf, "GrackleHeIII");
            break;
        case IO_grH2I:
            strcpy(buf, "GrackleH2I");
            break;
        case IO_grH2II:
            strcpy(buf, "GrackleH2II");
            break;
        case IO_grDI:
            strcpy(buf, "GrackleDI");
            break;
        case IO_grDII:
            strcpy(buf, "GrackleDII");
            break;
        case IO_grHDI:
            strcpy(buf, "GrackleHDI");
            break;
        case IO_fH2:
	    strcpy(buf, "FractionH2");
	        break;         
        case IO_NWINDLAUNCHES:
            strcpy(buf, "NWindLaunches");
            break;
        case IO_FOF_GROUPID:
            strcpy(buf, "HaloID");
            break;
        case IO_G0LOCAL:
            strcpy(buf, "LocalG0");
            break;
        case IO_TDUST:
            strcpy(buf, "DustTemperature");
            break;
        case IO_ACCKEY:
            strcpy(buf, "AccKey");
	        break;
        case IO_LASTENTRY:
            exit(2);
            break;
    }
}