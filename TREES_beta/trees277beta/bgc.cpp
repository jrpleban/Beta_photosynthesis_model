//***********************************************************
//Implementation of a BGC class
//
//This code is part of the TREES model
//  but a small effort needed to adapt to other models
// 2015, 2016 DSM
//
// This object manages states and fluxes associated with:
//   - allocation of carbon to structural pools
//   - allocation of non-structural carbon pools
//   - carbon and nitrogen cycles
//   - rhizosphere C and N dynamics
//
//******************** DSM 2015 *****************************

#include "simulator2.h"


BiogeochemicalCycles::BiogeochemicalCycles()
{
	leafResidueCarbon = NULL;
	leafResidueNitrogen = NULL;
	stemResidueCarbon = NULL;
	stemResidueNitrogen = NULL;
	rootResidueCarbon = NULL;
	rootResidueNitrogen = NULL;
	dCrootResidue = NULL;
	dNrootResidue = NULL;
	humusCarbon = NULL;
	humusNitrogen = NULL;
	soilAmmoniumNitrogen = NULL;
	soilNitrateNitrogen = NULL;
	rhizosphereCl = NULL;
	rhizosphereNl = NULL;
	rhizosphereLiveMicrobialCarbon = NULL;
	rhizosphereDeadMicrobialCarbon = NULL;
	rhizosphereLabileCarbon = NULL;
	rhizosphereMicrobialNitrogen = NULL;
	rhizosphereMineralNitrogen = NULL;
	rhizosphereAmmoniumNitrogen = NULL;
	rhizosphereNitrateNitrogen = NULL;
	fineRootBiomassCarbon = NULL;
	fineRootBiomassNitrogen = NULL;
	coarseRootBiomassCarbon = NULL;
	coarseRootBiomassNitrogen = NULL;
	rootNSC = NULL;
	rootMineralNitrogen = NULL;
	rootArea = NULL;
	liveStemCarbon = NULL;
	stemNSC = NULL;
	liveStemNitrogen = NULL;
	deadStemCarbon = NULL;
	deadStemNitrogen = NULL;
	leafBiomassCarbon = NULL;
	leafBiomassNitrogen = NULL;
	leafNSC = NULL;
	chloroplastStarch = NULL;
	chloroplastSugar = NULL;
	leafStoredNitrogen = NULL;
	leafRubiscoNitrogen = NULL;
	fruitCarbon = NULL;
	fruitNitrogen = NULL;
	plantNstatus = NULL;
	nitrogenLeaching = NULL;
	lat_Root_b_value_init = NULL;
	lat_Root_c_value_init = NULL;
	heterotrophicRespiration = NULL;
}

BiogeochemicalCycles::BiogeochemicalCycles(trees_params& treesParams)
{
//Number of soil-root layers not to exceed ULAT
	nRoots = treesParams.rmodules;
//Number of root sizes or orders
	nFineRootOrders = 5;
	nRootOrders = 10;

	allocatePools();
	initializePools(treesParams);
}

BiogeochemicalCycles::~BiogeochemicalCycles()
{
	clearPools();
}

//allocate memory to store BGC pools
void BiogeochemicalCycles::allocatePools()
{
	assert(nRoots > 0);
	assert(nRootOrders > 0);

//allocate pool for leaf residue C and N
	leafResidueCarbon = new double;
	assert(leafResidueCarbon != NULL);
	leafResidueNitrogen = new double;
	assert(leafResidueNitrogen != NULL);

//allocate pool for stem residue C and N
	stemResidueCarbon = new double;
	assert(stemResidueCarbon != NULL);
	stemResidueNitrogen = new double;
	assert(stemResidueNitrogen != NULL);

//allocate pools for root residue C and N, assuming we will lump fine and coarse root necromass
	rootResidueCarbon = new double*[nRoots];
	assert(rootResidueCarbon != NULL);
	rootResidueNitrogen = new double*[nRoots];
	assert(rootResidueNitrogen != NULL);
	dCrootResidue = new double*[nRoots];
	assert(dCrootResidue != NULL);
	dNrootResidue = new double*[nRoots];
	assert(dNrootResidue != NULL);

	for (int j = 0; j < nRoots; j++)
	{
		rootResidueCarbon[j] = new double[nRootOrders];
		assert(rootResidueCarbon[j] != NULL);
		rootResidueNitrogen[j] = new double[nRootOrders];
		assert(rootResidueNitrogen[j] != NULL);
		dCrootResidue[j] = new double[nRootOrders];
		assert(dCrootResidue[j] != NULL);
		dNrootResidue[j] = new double[nRootOrders];
		assert(dNrootResidue[j] != NULL);
	}

//allocate pool for recalcitrant (humus) C
	humusCarbon = new double[nRoots];
	assert(humusCarbon != NULL);
	humusNitrogen = new double[nRoots];
	assert(humusNitrogen != NULL);
	soilAmmoniumNitrogen = new double[nRoots];
	assert(soilAmmoniumNitrogen != NULL);
	soilNitrateNitrogen = new double[nRoots];
	assert(soilNitrateNitrogen != NULL);

//State variables for litter C and litter N
	rhizosphereCl = new double*[nRoots];
	assert(rhizosphereCl != NULL);
	rhizosphereNl = new double*[nRoots];
	assert(rhizosphereNl != NULL);

//allocate microbial pools
	rhizosphereLiveMicrobialCarbon = new double*[nRoots];
	assert(rhizosphereLiveMicrobialCarbon != NULL);
	rhizosphereDeadMicrobialCarbon = new double*[nRoots];
	assert(rhizosphereDeadMicrobialCarbon != NULL);
	rhizosphereLabileCarbon = new double*[nRoots];
	assert(rhizosphereLabileCarbon != NULL);
	rhizosphereMicrobialNitrogen = new double*[nRoots];
	assert(rhizosphereMicrobialNitrogen != NULL);
	rhizosphereMineralNitrogen = new double*[nRoots];
	assert(rhizosphereMineralNitrogen != NULL);
	rhizosphereAmmoniumNitrogen = new double*[nRoots];
	assert(rhizosphereMineralNitrogen != NULL);
	rhizosphereNitrateNitrogen = new double*[nRoots];
	assert(rhizosphereMineralNitrogen != NULL);

	for (int j = 0; j < nRoots; j++)
	{
		rhizosphereCl[j] = new double[nRootOrders];
		rhizosphereNl[j] = new double[nRootOrders];
		rhizosphereLiveMicrobialCarbon[j] = new double[nRootOrders];
		assert(rhizosphereLiveMicrobialCarbon[j] != NULL);
		rhizosphereDeadMicrobialCarbon[j] = new double[nRootOrders];
		assert(rhizosphereDeadMicrobialCarbon[j] != NULL);
		rhizosphereLabileCarbon[j] = new double[nRootOrders];
		assert(rhizosphereLabileCarbon[j] != NULL);
		rhizosphereMicrobialNitrogen[j] = new double[nRootOrders];
		assert(rhizosphereMicrobialNitrogen[j] != NULL);
		rhizosphereMineralNitrogen[j] = new double[nRootOrders];
		assert(rhizosphereMineralNitrogen[j] != NULL);
		rhizosphereAmmoniumNitrogen[j] = new double[nRootOrders];
		assert(rhizosphereAmmoniumNitrogen[j] != NULL);
		rhizosphereNitrateNitrogen[j] = new double[nRootOrders];
		assert(rhizosphereNitrateNitrogen[j] != NULL);
	}

//allocate heterotrophic respiration
	heterotrophicRespiration = new double[nRoots];
	assert(heterotrophicRespiration != NULL);

//allocate root pools
	fineRootBiomassCarbon = new double*[nRoots];
	assert(fineRootBiomassCarbon != NULL);
	fineRootBiomassNitrogen = new double*[nRoots];
	assert(fineRootBiomassNitrogen != NULL);
	coarseRootBiomassCarbon = new double*[nRoots];
	assert(coarseRootBiomassCarbon != NULL);
	coarseRootBiomassNitrogen = new double*[nRoots];
	assert(coarseRootBiomassNitrogen != NULL);
	rootNSC = new double*[nRoots];
	assert(rootNSC != NULL);
	rootMineralNitrogen = new double*[nRoots];
	assert(rootMineralNitrogen != NULL);
	rootArea = new double*[nRoots];
	assert(rootArea != NULL);

	for (int j = 0; j < nRoots; j++)
	{
		fineRootBiomassCarbon[j] = new double[nRootOrders];
		assert(fineRootBiomassCarbon[j] != NULL);
		fineRootBiomassNitrogen[j] = new double[nRootOrders];
		assert(fineRootBiomassNitrogen[j] != NULL);
		coarseRootBiomassCarbon[j] = new double[nRootOrders];
		assert(coarseRootBiomassCarbon[j] != NULL);
		coarseRootBiomassNitrogen[j] = new double[nRootOrders];
		assert(coarseRootBiomassNitrogen[j] != NULL);
		rootNSC[j] = new double[nRootOrders];
		assert(rootNSC[j] != NULL);
		rootMineralNitrogen[j] = new double[nRootOrders];
		assert(rootMineralNitrogen[j] != NULL);
		rootArea[j] = new double[nRootOrders];
		assert(rootArea[j] != NULL);
	}

//allocate stem pools
	liveStemCarbon = new double;
	assert(liveStemCarbon != NULL);
	stemNSC = new double;
	assert(stemNSC != NULL);
	liveStemNitrogen = new double;
	assert(liveStemNitrogen != NULL);
	deadStemCarbon = new double;
	assert(deadStemCarbon != NULL);
	deadStemNitrogen = new double;
	assert(deadStemNitrogen != NULL);

//allocate leaf pools
	leafBiomassCarbon = new double;
	assert(leafBiomassCarbon != NULL);
	leafBiomassNitrogen = new double;
	assert(leafBiomassNitrogen != NULL);
	leafNSC = new double;
	assert(leafNSC != NULL);
	chloroplastStarch = new double;
	assert(chloroplastStarch != NULL);
	chloroplastSugar = new double;
	assert(chloroplastSugar != NULL);
	leafStoredNitrogen = new double;
	assert(leafStoredNitrogen != NULL);
	leafRubiscoNitrogen = new double;
	assert(leafRubiscoNitrogen != NULL);

//allocation reproductive pools
	fruitCarbon = new double;
	assert(fruitCarbon != NULL);
	fruitNitrogen = new double;
	assert(fruitNitrogen != NULL);

//allocation for plant nitrogen status
	plantNstatus = new double;
	assert(plantNstatus != NULL);

//allocation for the nitrogen leaching flux variable
	nitrogenLeaching = new double[nRoots];
	assert(nitrogenLeaching != NULL);

//allocation for lateral root Weibulls
	lat_Root_b_value_init = new double;
	assert(lat_Root_b_value_init != NULL);
	lat_Root_c_value_init = new double;
	assert(lat_Root_c_value_init != NULL);
}

//set BGC pools to initial values
//the assumption with this code is that many state variables are correlated
//   and can be estimates from a few input paramaters such as SLA and leaf lifespan
//Common C:N for meterials:
//    > deciduous leaves 32
//        > litter       64
//    > needle leaf     110
//        > litter      220
//    > microbes          8
//    > soil humus       10
void BiogeochemicalCycles::initializePools(trees_params treesParams)
{
	double totalRootArea, scalar, rootScalar, sumCl, leafN;
	double SLAscalar, SLAreference, totDepth;
	double CN;

	SLAreference = 22.0;
	SLAscalar = SLAreference / treesParams.SLA;
	if (SLAscalar > 3.0)
	{
		SLAscalar = 3.0;
	}
	if (SLAscalar < 1.0)
	{
		SLAscalar = 1.0;
	}

	totDepth = 0.0;
	for (int j = 0; j < treesParams.rmodules; j++)
	{
		totDepth += treesParams.drax[j+3];
	}

//assume initial leaf residue is proportion to leaf biomass and increasing with declining SLA
//units are kgC ha-1 and kgN ha-1
//set leaf residue N to 50% of live leaf N, assuming retranslocation
	if (treesParams.usePhenology == true) //best for perennial plants
	{
        	leafResidueCarbon[0] = (1.0/treesParams.leafLifeSpan)*treesParams.lai / treesParams.SLA * 10000.0 * SLAscalar;
        	leafResidueNitrogen[0] = (0.5 * treesParams.Nleaf * treesParams.lai * 10000.0)*SLAscalar;
//set stem residue to 5% of standing stem biomass
//stem [N] is assumed to be low (0.3% of carbon)
        	stemResidueCarbon[0] = 0.05 * treesParams.Cstem;
        	stemResidueNitrogen[0] = 0.003 * stemResidueCarbon[0];
	}
	else //best for annual plants
	{
		leafResidueCarbon[0] = 0.00001;
		leafResidueNitrogen[0] = 0.00000001;
        	stemResidueCarbon[0] = 0.00001;
        	stemResidueNitrogen[0] = 0.00000001;
	}


//set root and rhizosphere state variables
//   there are nRoots = number of soil-root layers in TREES, defined in "param_mod" input file
//   there are 10 root orders, assuming the first 5 are fine root (from 1/4 mm up to 4 mm diameter)
        for (int j = 0; j < nRoots; j++)
        {
		scalar = 8.0;
		rootScalar = 1.0;
		sumCl = 0.0;
		totDepth += treesParams.drax[j+3];

		heterotrophicRespiration[j] = 0.0;

		for (int k = 0; k < nFineRootOrders; k++)
		{
//State variables for litter C and litter N
                	fineRootBiomassCarbon[j][k] = treesParams.Croot * treesParams.ar[j+3] * 0.2;
                	fineRootBiomassNitrogen[j][k] = fineRootBiomassCarbon[j][k] / (20.0*rootScalar*SLAscalar);
                	coarseRootBiomassCarbon[j][k] = 0.0;
                	coarseRootBiomassNitrogen[j][k] = coarseRootBiomassCarbon[j][k] / (60.0*SLAscalar);

                	rhizosphereCl[j][k] = (scalar * SLAscalar)*
				(fineRootBiomassCarbon[j][k] + 0.1*coarseRootBiomassCarbon[j][k]);
			sumCl += rhizosphereCl[j][k];
                	rhizosphereNl[j][k] = rhizosphereCl[j][k] / (20.0+20.0*(1.0-1.0/SLAscalar));
			rootResidueCarbon[j][k] = scalar*fineRootBiomassCarbon[j][k];
			rootResidueNitrogen[j][k] = scalar*fineRootBiomassNitrogen[j][k];
                	rhizosphereLiveMicrobialCarbon[j][k] = 0.06*rhizosphereCl[j][k];
                	rhizosphereDeadMicrobialCarbon[j][k] = rhizosphereLiveMicrobialCarbon[j][k];
			sumCl += rhizosphereLiveMicrobialCarbon[j][k] + rhizosphereDeadMicrobialCarbon[j][k];
                	rhizosphereMicrobialNitrogen[j][k] = rhizosphereLiveMicrobialCarbon[j][k] / (8.0+2.0*(1.0-1.0/SLAscalar));
                	rhizosphereMineralNitrogen[j][k] = 1.0*rhizosphereDeadMicrobialCarbon[j][k] / (8.0+2.0*(1.0-1.0/SLAscalar));
                	rhizosphereAmmoniumNitrogen[j][k] = 0.1 * rhizosphereMineralNitrogen[j][k];
                	rhizosphereNitrateNitrogen[j][k] = 0.9 * rhizosphereMineralNitrogen[j][k];

                	rootNSC[j][k] = 0.5*treesParams.leafNSCscalar * fineRootBiomassCarbon[j][k];
                	rhizosphereLabileCarbon[j][k] = 0.1*rootNSC[j][k];
			CN = fineRootBiomassCarbon[j][k]/fineRootBiomassNitrogen[j][k];
			//rootMineralNitrogen[j][k] = 2.0*rootNSC[j][k]/CN;
			rootMineralNitrogen[j][k] = fineRootBiomassNitrogen[j][k];
                	//rootMineralNitrogen[j][k] = 1.5/SLAscalar * fineRootBiomassNitrogen[j][k] * (1.0/treesParams.leafLifeSpan);;
//allow for a scaling of the initial nutrient levels in the microbiome
                        if (treesParams.microbiomeScalar < 0.0)
                        {
				treesParams.microbiomeScalar = 0.0;
			}
                        rhizosphereMineralNitrogen[j][k] *= treesParams.microbiomeScalar;
			rhizosphereAmmoniumNitrogen[j][k] *= treesParams.microbiomeScalar;
			rhizosphereNitrateNitrogen[j][k] *= treesParams.microbiomeScalar;
//                        rootMineralNitrogen[j][k] *= treesParams.microbiomeScalar;
			rhizosphereLiveMicrobialCarbon[j][k] *= treesParams.microbiomeScalar;
			rhizosphereDeadMicrobialCarbon[j][k] *= treesParams.microbiomeScalar;
			rhizosphereMicrobialNitrogen[j][k] *= treesParams.microbiomeScalar;

			scalar *= 0.50;
			rootScalar *= 1.25;
		}
		for (int k = nFineRootOrders; k < nRootOrders; k++)
		{

                	fineRootBiomassCarbon[j][k] = 0.0;
                	fineRootBiomassNitrogen[j][k] = fineRootBiomassCarbon[j][k] / (30.0*SLAscalar);
                	coarseRootBiomassCarbon[j][k] = treesParams.Croot_coarse * treesParams.ar[j+3] * 0.16;
                	coarseRootBiomassCarbon[j][k] += treesParams.Croot_coarse * treesParams.drax[j+3]/ totDepth * 0.04;
                	coarseRootBiomassNitrogen[j][k] = coarseRootBiomassCarbon[j][k] / (60.0*SLAscalar);

                	rhizosphereCl[j][k] = (scalar * SLAscalar)*
					(fineRootBiomassCarbon[j][k] + coarseRootBiomassCarbon[j][k]);
			sumCl += rhizosphereCl[j][k];
                	rhizosphereNl[j][k] = rhizosphereCl[j][k] / (20.0+20.0*(1.0-1.0/SLAscalar));
			rootResidueCarbon[j][k] = scalar*coarseRootBiomassCarbon[j][k];
			rootResidueNitrogen[j][k] = scalar*coarseRootBiomassNitrogen[j][k];
                	rhizosphereLiveMicrobialCarbon[j][k] = 0.06 * rhizosphereCl[j][k];
                	rhizosphereDeadMicrobialCarbon[j][k] = rhizosphereLiveMicrobialCarbon[j][k];
			sumCl += rhizosphereLiveMicrobialCarbon[j][k] + rhizosphereDeadMicrobialCarbon[j][k];
                	rhizosphereMicrobialNitrogen[j][k] = rhizosphereLiveMicrobialCarbon[j][k] / (8.0+2.0*(1.0-1.0/SLAscalar));
                	rhizosphereMineralNitrogen[j][k] = 1.0*rhizosphereDeadMicrobialCarbon[j][k] / (8.0+2.0*(1.0-1.0/SLAscalar));
                	rhizosphereAmmoniumNitrogen[j][k] = 0.1 * rhizosphereMineralNitrogen[j][k];
                	rhizosphereNitrateNitrogen[j][k] = 0.9 * rhizosphereMineralNitrogen[j][k];

                	rootNSC[j][k] = 0.5*treesParams.leafNSCscalar * coarseRootBiomassCarbon[j][k];
                	rhizosphereLabileCarbon[j][k] = 0.1*rootNSC[j][k];
			CN = coarseRootBiomassCarbon[j][k]/coarseRootBiomassNitrogen[j][k];
			//rootMineralNitrogen[j][k] = 2.0*rootNSC[j][k]/CN;
			rootMineralNitrogen[j][k] = fineRootBiomassNitrogen[j][k];
                	//rootMineralNitrogen[j][k] = 1.5/SLAscalar * coarseRootBiomassNitrogen[j][k] * (1.0/treesParams.leafLifeSpan);;
//allow for a scaling of the initial nutrient levels in the microbiome
                        if (treesParams.microbiomeScalar < 0.0)
                        {
				treesParams.microbiomeScalar = 0.0;
			}
                        rhizosphereMineralNitrogen[j][k] *= treesParams.microbiomeScalar;
			rhizosphereAmmoniumNitrogen[j][k] *= treesParams.microbiomeScalar;
			rhizosphereNitrateNitrogen[j][k] *= treesParams.microbiomeScalar;
//                        rootMineralNitrogen[j][k] *= treesParams.microbiomeScalar;
			rhizosphereLiveMicrobialCarbon[j][k] *= treesParams.microbiomeScalar;
			rhizosphereDeadMicrobialCarbon[j][k] *= treesParams.microbiomeScalar;
			rhizosphereMicrobialNitrogen[j][k] *= treesParams.microbiomeScalar;

			scalar *= 0.50;
		}
//State variables for stabilized soil carbon
                humusCarbon[j] = treesParams.Csoil * treesParams.ar[j+3] - sumCl;
                humusNitrogen[j] = 0.09 * humusCarbon[j];

//State variables for mineral nitrogen in the bulk soil
		soilAmmoniumNitrogen[j] = 0.1*humusNitrogen[j];
		soilNitrateNitrogen[j] = 0.9*humusNitrogen[j];
//Flux variable for nitrogen leaching set to zero
		nitrogenLeaching[j] = 0.0;

        }

	totalRootArea = computeRootArea(treesParams);

        liveStemCarbon[0] = treesParams.Csapwood;
        stemNSC[0] = 0.5*treesParams.leafNSCscalar*0.4 * liveStemCarbon[0];
        liveStemNitrogen[0] = 0.003 * liveStemCarbon[0];
        deadStemCarbon[0] = treesParams.Cstem;
        deadStemNitrogen[0] = 0.003 * deadStemCarbon[0];

	if (treesParams.leafLifeSpan > 1.0)
	{
        	leafBiomassCarbon[0] = (1.0 - 1.0/treesParams.leafLifeSpan)*treesParams.lai/treesParams.SLA*10000.0+0.01;
	}
	else
	{
        	leafBiomassCarbon[0] = (1.0 - 0.9/treesParams.leafLifeSpan)*treesParams.lai/treesParams.SLA*10000.0+0.01;
	}
        leafNSC[0] = 0.5*leafBiomassCarbon[0]*treesParams.leafNSCscalar;
	leafN = treesParams.Nleaf*treesParams.lai*10000.0;
        leafBiomassNitrogen[0] = max(0.1,(treesParams.leafLifeSpan-1.0)/treesParams.leafLifeSpan)*leafN;
	leafRubiscoNitrogen[0] = leafBiomassNitrogen[0] * treesParams.Nrubisco;
	chloroplastStarch[0] = 0.01*leafRubiscoNitrogen[0]*leafBiomassCarbon[0]/leafBiomassNitrogen[0];
	chloroplastSugar[0] = 0.0;
	CN = leafBiomassCarbon[0]/leafBiomassNitrogen[0];
	//leafStoredNitrogen[0] = 2.0*(leafNSC[0]+stemNSC[0])/CN;
	leafStoredNitrogen[0] = 0.5*leafBiomassNitrogen[0];
	if (treesParams.lai < 1.0 && treesParams.lai > 0.0)
	{
		leafStoredNitrogen[0] /= treesParams.lai;
	}
	//leafStoredNitrogen[0] = 5.0/treesParams.leafLifeSpan*leafBiomassNitrogen[0];

        //leafStoredNitrogen[0] *= treesParams.microbiomeScalar;

        fruitCarbon[0] = 0.0;
        fruitNitrogen[0] = 0.0;

	plantNstatus[0] = 1.0; //values are between 0 and 1
//keep a record of the input lateral root Weibull curve parameters
	lat_Root_b_value_init[0] = treesParams.lat_Root_b_value;
	lat_Root_c_value_init[0] = treesParams.lat_Root_c_value;
}


//clear and delete BGC pool variables
void BiogeochemicalCycles::clearPools()
{
	if (leafResidueCarbon != NULL)
	{
		delete leafResidueCarbon;
		leafResidueCarbon = NULL;
	}
	if (leafResidueNitrogen != NULL)
	{
		delete leafResidueNitrogen;
		leafResidueNitrogen = NULL;
	}
	if (stemResidueCarbon != NULL)
	{
		delete stemResidueCarbon;
		stemResidueCarbon = NULL;
	}
	if (stemResidueNitrogen != NULL)
	{
		delete stemResidueNitrogen;
		stemResidueNitrogen = NULL;
	}
	if (rootResidueCarbon != NULL)
	{
		for (int j = 0; j < nRoots; j++)
		{
			if (rootResidueCarbon[j] != NULL)
			{
				delete[] rootResidueCarbon[j];
				rootResidueCarbon[j] = NULL;
			}
		}
		delete[] rootResidueCarbon;
		rootResidueCarbon = NULL;
	}
	if (rootResidueNitrogen != NULL)
	{
		for (int j = 0; j < nRoots; j++)
		{
			if (rootResidueNitrogen[j] != NULL)
			{
				delete[] rootResidueNitrogen[j];
				rootResidueNitrogen[j] = NULL;
			}
		}
		delete[] rootResidueNitrogen;
		rootResidueNitrogen = NULL;
	}
	if (humusCarbon != NULL)
	{
		delete[] humusCarbon;
		humusCarbon = NULL;
	}
	if (humusNitrogen != NULL)
	{
		delete[] humusNitrogen;
		humusNitrogen = NULL;
	}
	if (heterotrophicRespiration != NULL)
	{
		delete[] heterotrophicRespiration;
		heterotrophicRespiration = NULL;
	}
	if (rhizosphereCl != NULL)
	{
		for (int j = 0; j < nRoots; j++)
		{
			if (rhizosphereCl[j] != NULL)
			{
				delete[] rhizosphereCl[j];
				rhizosphereCl[j] = NULL;
			}
		}
		delete[] rhizosphereCl;
		rhizosphereCl = NULL;
	}
	if (rhizosphereNl != NULL)
	{
		for (int j = 0; j < nRoots; j++)
		{
			if (rhizosphereNl[j] != NULL)
			{
				delete[] rhizosphereNl[j];
				rhizosphereNl[j] = NULL;
			}
		}
		delete[] rhizosphereNl;
		rhizosphereNl = NULL;
	}
	if (rhizosphereLiveMicrobialCarbon != NULL)
	{
		for (int j = 0; j < nRoots; j++)
		{
			if (rhizosphereLiveMicrobialCarbon[j] != NULL)
			{
				delete[] rhizosphereLiveMicrobialCarbon[j];
				rhizosphereLiveMicrobialCarbon[j] = NULL;
			}
		}
		delete[] rhizosphereLiveMicrobialCarbon;
		rhizosphereLiveMicrobialCarbon = NULL;
	}
	if (rhizosphereDeadMicrobialCarbon != NULL)
	{
		for (int j = 0; j < nRoots; j++)
		{
			if (rhizosphereDeadMicrobialCarbon[j] != NULL)
			{
				delete[] rhizosphereDeadMicrobialCarbon[j];
				rhizosphereDeadMicrobialCarbon[j] = NULL;
			}
		}
		delete[] rhizosphereDeadMicrobialCarbon;
		rhizosphereDeadMicrobialCarbon = NULL;
	}
	if (rhizosphereLabileCarbon != NULL)
	{
		for (int j = 0; j < nRoots; j++)
		{
			if (rhizosphereLabileCarbon[j] != NULL)
			{
				delete[] rhizosphereLabileCarbon[j];
				rhizosphereLabileCarbon[j] = NULL;
			}
		}
		delete[] rhizosphereLabileCarbon;
		rhizosphereLabileCarbon = NULL;
	}
	if (rhizosphereMicrobialNitrogen != NULL)
	{
		for (int j = 0; j < nRoots; j++)
		{
			if (rhizosphereMicrobialNitrogen[j] != NULL)
			{
				delete[] rhizosphereMicrobialNitrogen[j];
				rhizosphereMicrobialNitrogen[j] = NULL;
			}
		}
		delete[] rhizosphereMicrobialNitrogen;
		rhizosphereMicrobialNitrogen = NULL;
	}
	if (rhizosphereMineralNitrogen != NULL)
	{
		for (int j = 0; j < nRoots; j++)
		{
			if (rhizosphereMineralNitrogen[j] != NULL)
			{
				delete[] rhizosphereMineralNitrogen[j];
				rhizosphereMineralNitrogen[j] = NULL;
			}
		}
		delete[] rhizosphereMineralNitrogen;
		rhizosphereMineralNitrogen = NULL;
	}
	if (rhizosphereAmmoniumNitrogen != NULL)
	{
		for (int j = 0; j < nRoots; j++)
		{
			if (rhizosphereAmmoniumNitrogen[j] != NULL)
			{
				delete[] rhizosphereAmmoniumNitrogen[j];
				rhizosphereAmmoniumNitrogen[j] = NULL;
			}
		}
		delete[] rhizosphereAmmoniumNitrogen;
		rhizosphereAmmoniumNitrogen = NULL;
	}
	if (rhizosphereNitrateNitrogen != NULL)
	{
		for (int j = 0; j < nRoots; j++)
		{
			if (rhizosphereNitrateNitrogen[j] != NULL)
			{
				delete[] rhizosphereNitrateNitrogen[j];
				rhizosphereNitrateNitrogen[j] = NULL;
			}
		}
		delete[] rhizosphereNitrateNitrogen;
		rhizosphereNitrateNitrogen = NULL;
	}
	if (fineRootBiomassCarbon != NULL)
	{
		for (int j = 0; j < nRoots; j++)
		{
			if (fineRootBiomassCarbon[j] != NULL)
			{
				delete[] fineRootBiomassCarbon[j];
				fineRootBiomassCarbon[j] = NULL;
			}
		}
		delete[] fineRootBiomassCarbon;
		fineRootBiomassCarbon = NULL;
	}
	if (fineRootBiomassNitrogen != NULL)
	{
		for (int j = 0; j < nRoots; j++)
		{
			if (fineRootBiomassNitrogen[j] != NULL)
			{
				delete[] fineRootBiomassNitrogen[j];
				fineRootBiomassNitrogen[j] = NULL;
			}
		}
		delete[] fineRootBiomassNitrogen;
		fineRootBiomassNitrogen = NULL;
	}
	if (coarseRootBiomassCarbon != NULL)
	{
		for (int j = 0; j < nRoots; j++)
		{
			if (coarseRootBiomassCarbon[j] != NULL)
			{
				delete[] coarseRootBiomassCarbon[j];
				coarseRootBiomassCarbon[j] = NULL;
			}
		}
		delete[] coarseRootBiomassCarbon;
		coarseRootBiomassCarbon = NULL;
	}
	if (coarseRootBiomassNitrogen != NULL)
	{
		for (int j = 0; j < nRoots; j++)
		{
			if (coarseRootBiomassNitrogen[j] != NULL)
			{
				delete[] coarseRootBiomassNitrogen[j];
				coarseRootBiomassNitrogen[j] = NULL;
			}
		}
		delete[] coarseRootBiomassNitrogen;
		coarseRootBiomassNitrogen = NULL;
	}
	if (rootNSC != NULL)
	{
		for (int j = 0; j < nRoots; j++)
		{
			if (rootNSC[j] != NULL)
			{
				delete[] rootNSC[j];
				rootNSC[j] = NULL;
			}
		}
		delete[] rootNSC;
		rootNSC = NULL;
	}
	if (rootMineralNitrogen != NULL)
	{
		for (int j = 0; j < nRoots; j++)
		{
			if (rootMineralNitrogen[j] != NULL)
			{
				delete[] rootMineralNitrogen[j];
				rootMineralNitrogen[j] = NULL;
			}
		}
		delete[] rootMineralNitrogen;
		rootMineralNitrogen = NULL;
	}
	if (rootArea != NULL)
	{
		for (int j = 0; j < nRoots; j++)
		{
			if (rootArea[j] != NULL)
			{
				delete[] rootArea[j];
				rootArea[j] = NULL;
			}
		}
		delete[] rootArea;
		rootArea = NULL;
	}
	if (liveStemCarbon != NULL)
	{
		delete liveStemCarbon;
		liveStemCarbon = NULL;
	}
	if (stemNSC != NULL)
	{
		delete stemNSC;
		stemNSC = NULL;
	}
	if (liveStemNitrogen != NULL)
	{
		delete liveStemNitrogen;
		liveStemNitrogen = NULL;
	}
	if (deadStemCarbon != NULL)
	{
		delete deadStemCarbon;
		deadStemCarbon = NULL;
	}
	if (deadStemNitrogen != NULL)
	{
		delete deadStemNitrogen;
		deadStemNitrogen = NULL;
	}
	if (leafBiomassCarbon != NULL)
	{
		delete leafBiomassCarbon;
		leafBiomassCarbon = NULL;
	}
	if (leafBiomassNitrogen != NULL)
	{
		delete leafBiomassNitrogen;
		leafBiomassNitrogen = NULL;
	}
	if (leafNSC != NULL)
	{
		delete leafNSC;
		leafNSC = NULL;
	}
	if (chloroplastStarch != NULL)
	{
		delete chloroplastStarch;
		chloroplastStarch = NULL;
	}
	if (chloroplastSugar != NULL)
	{
		delete chloroplastSugar;
		chloroplastSugar = NULL;
	}
	if (leafStoredNitrogen != NULL)
	{
		delete leafStoredNitrogen;
		leafStoredNitrogen = NULL;
	}
	if (leafRubiscoNitrogen != NULL)
	{
		delete leafRubiscoNitrogen;
		leafRubiscoNitrogen = NULL;
	}
	if (fruitCarbon != NULL)
	{
		delete fruitCarbon;
		fruitCarbon = NULL;
	}
	if (fruitNitrogen != NULL)
	{
		delete fruitNitrogen;
		fruitNitrogen = NULL;
	}
	if (plantNstatus != NULL)
	{
		delete plantNstatus;
		plantNstatus = NULL;
	}
	if (nitrogenLeaching != NULL)
	{
		delete nitrogenLeaching;
		nitrogenLeaching = NULL;
	}
	if (lat_Root_b_value_init != NULL)
	{
		delete lat_Root_b_value_init;
		lat_Root_b_value_init = NULL;
	}
	if (lat_Root_c_value_init != NULL)
	{
		delete lat_Root_c_value_init;
		lat_Root_c_value_init = NULL;
	}
}


//
// photosyntheis()
//     This function computes leaf photosynthesis for C3 plants,
//	returning also stomatal conductance
//
//     August 2017 DSM - Modification to allow noctural stomatal conductance
//
double BiogeochemicalCycles::photosynthesis(trees_params treesParams,
                                            double Jmax_mult,
                                            double thetaJ,
                                            double phiJ,
                                            double daily_md_lwp,
                                            struct farqin in,
                                            struct farqout& out)
{
	double t;      			// (deg C) temperature
	double tk;     			// (K) absolute temperature
	double g, gin;      		// (umol/m2/s/Pa) conductance to CO2
	double O2;     			// (Pa) atmospheric partial pressure O2
	double Ca;     			// (Pa) atmospheric partial pressure CO2
	double gammaStar;  		// (Pa) co2 compensation point, no dark respiration
	double Kc25, Kc;     		// (Pa) MM constant for carboxylase reaction
	double q10Kc;
	double Ko25, Ko;     		// (Pa) MM constant for oxygenase reaction
	double q10Ko;
	double act25, act;    		// (umol/kgRubisco/s) Rubisco activity
	double q10act;
	double Rd;     			// (umol/m2/s) dark respiration rate
	double Vcmax, Vcmax25;   	// (umol/m2/s) maximum carboxylation velocity
	double Jmax, Jmax25;   		// (umol/m2/s) maximum rate of electron transport
	double J;      			// (umol/m2/s) maximum rate of Rubisco regeneration
	double Av;     			// (umol/m2/s) Rubisco limited assimilation rate
	double Aj;     			// (umol/m2/s) RuBP regeneration limited assim rate
	double As;      		// (umol/m2/s) Sink-limited assimilation rate
	double A;      			// (umol/m2/s) net assimilation rate
	double A_kg;			// (kgC ha-1)
	//double Ile;   		// (umol/m2/s) PAR effectively absorbed by PSII per unit leaf area
	double Ji;      		// ((umol/m2/s) rate of electron transport per unit absorbed PAR
    double phi2;
    double phi2_300;
    double betaA;
    double alpha;
    double s;        // s is lumper parameter to estimate PAR abs by PSII (abs*alt e- flow * PSII/PSI partitioning)
    double kappa;    // minimum phi2 at saturating light
    // s and kappa are currently fixed based on population level median posterior distributinos of betaA decay photosynthesis model
    double m_phi;    // slope of phi2 against increasly (-) leaf psi
    // m_phi is comes from R500 mortality data set (Guadagno, PlantPhys 2017)
    double phi2_hyd;
	double aa,bb,cc,det;
	double kappa1, kappa2, kappa3;	//constants used in Katul et al, 2003, Plant, Cell and Env
	double alpha1, alpha2;   	//constants used in Katul et al, 2003, Plant, Cell and Env
	double Ci, D;
	double E;       		//transpiration used to check validity of Ci
	double g_mult;
	double Rd_mult = treesParams.Rd_mult;
    double phiJ_shd = treesParams.phiJ_shd;
    double phiJ_sun = treesParams.phiJ_sun;
/* local variables  */

        t = in.t;               	//temperature should be leaf temperature
        tk = t + 273.15;        	//leaf temperature (deg. K)
        g = gin = in.g*1.0e6/(R*tk);    //convert conductance from m/s -. umol/m2/s/Pa
        if ( g < 0.00000001 )
        {
                g = 0.00000001;
        }
        Ca = in.co2 * 1.0e-6 * in.pa;   //convert atmospheric CO2 from ppm to Pa
        O2 = 0.2095 * in.pa; 		//atmospheric O2 in Pa, assumes 20.95% O2 by volume

// correct kinetic constants for temperature, and do unit conversions
//Future work - replace with Arrhenius kinetics
        Kc25 = treesParams.Kc25;
        q10Kc = treesParams.q10Kc;
        Ko25 = treesParams.Ko25;
        q10Ko = treesParams.q10Ko;
        act25 = treesParams.act25;
        q10act = treesParams.q10act;
        Ko = Ko25 * pow(q10Ko, (t-25.0)/10.0);
        if (t > 15.0)
        {
                Kc = Kc25 * pow(q10Kc, (t-25.0)/10.0);
                act = act25 * pow(q10act, (t-25.0)/10.0);
        }
        else
        {
                Kc = Kc25 * pow(1.8*q10Kc, (t-15.0)/10.0) / q10Kc;
                act = act25 * pow(1.8*q10act, (t-15.0)/10.0) / q10act;
        }

        act = act * 1.0e6 / 60.0;     // umol/mg/min to umol/kg/s

// calculate gammaStar (Pa) - DePury and Farquhar, 1997
        gammaStar = 3.69+0.188*(t-25.0)+0.0036*(t-25.0)*(t-25.0);

// calculate Vcmax from leaf nitrogen data and Rubisco activity

/* kg Nleaf   kg NRub    kg Rub      umol            umol
           -------- X -------  X ------- X ---------   =   --------
              m2      kg Nleaf   kg NRub   kg Rub * s       m2 * s

             (lnc)  X  (flnr)  X  (fnr)  X   (act)     =    (Vcmax)
*/

        Vcmax = in.lnc * in.flnr * fnr * act;

/* Leaf respiration not including photorespiration calc. as
           Rd = 0.0089Vcmax
           from Leuning et al.  1995     */

        Rd = Rd_mult * Vcmax;
        out.Rd = Rd;

/* calculate Jmax = f(Vcmax), reference:
        Wullschleger, S.D., 1993.  Biochemical limitations to carbon assimilation
                in C3 plants - A retrospective analysis of the A/Ci curves from
                109 species. Journal of Experimental Botany, 44:907:920.
*/

//calculate J = f(Jmax, ppfd), reference: de Pury and Farquhar 1997 Plant Cell and Env.
//Jmax proportion to Vcmax at the reference temperature, 25 deg. C (JUNE 2008 DSM)
        Vcmax25 = in.lnc * in.flnr * fnr * act25 * 1.0e6/60.0;
        //Jmax25 = Jmax_mult*Vcmax25;

        double Ea = 37000.0;   //Ativation energy, kJ mol-1
        double S = 710.0;      //Electrong transport temperature response parameter, J K-1 mol-1
        double H = 220000.0;   //Electron transport temperature curvature parameter, J mol-1


    
       // Jmax = Jmax25*exp((tk-298.0)*Ea/(R*tk*298.0))*
		//		(1.0+exp((S*298.0-H)/(R*298.0)))/(1.0+exp((S*tk-H)/(R*tk)));

//phiJ is effective quantum yield, mol electrons mol-1 photons
//assumes phiJ = quantum yield (mol C mol-1 photons) * 4 e- per C * leaf absorptance
//phiJ can be varied between sun and shade leaves to take into consideration direct versus diffuse rad'n
//therefore, Ile no longer used, but embedded in phiJ
        if (in.irad > 2200.0)
        {
                in.irad = 2200.0;
        }
        if (in.irad < 0.0)
        {
                in.irad = 0.0;
        }
    /*
        Ji = in.irad * phiJ;

        aa = thetaJ;
        bb = -Ji -Jmax;
        cc = Ji*Jmax;
        J = (-bb - sqrt(bb*bb - 4.0*aa*cc))/(2.0*aa);
     */
    // modification to photosynthesis routine to describe estimate eletron rate using paratmers from chl fluoresnce and the realtionship between phi2 and leaf water potential (JRP 1/2018)
    // phi2_300 is PSII efficiecy under optimal conditions
    // 0.07 slope of phi2 vs leafpsi relationship based on R500 mortality data set (Guadagno, PlantPhys 2017)
    // kappa is minimum phi2 at saturating light
    // betaA is decay rate in phi2 underin increasing light
    // alpha is maximum phi2 efficienct (ligh adapted)
    // leafpsi_md is previous daya mid-day leafpsi from hydrualic model
    // s is lumper parameter to estimate PAR abs by PSII (abs*alt e- flow * PSII/PSI partitioning)
    //double md_lwp;
    //md_lwp= leafpsi;
    
    m_phi = Jmax_mult;  /// change par names here for simplicity from original trees  .p file
    phi2_300 = thetaJ;  /// change par names here for simplicity from original trees  .p file
    alpha = phiJ_sun;   /// change par names here for simplicity from original trees  .p file
    kappa = phiJ_shd; /// change par names here for simplicity from original trees  .p file
    
    /// beta decay implimentaion - version 1 beta decrease with leaf psi only
    
    phi2_hyd = m_phi*(daily_md_lwp) + phi2_300;
    betaA =(log((phi2_hyd-kappa)/(alpha-kappa)))/300;
    s= 0.31;
    //(betaA  * 45)+ 0.31;
    phi2 = (alpha-kappa)* (exp(betaA * (in.irad )))+kappa;
    J = in.irad * phi2* s * 2;
    
    //cout << "daily_md_lwp = " << daily_md_lwp << endl;
    //cout << "phi2_hyd = " << phi2_hyd << endl;
    //cout << "betaA = " << betaA << endl;
    //cout << "s = " << s << endl;

/* solve for Av and Aj using the quadratic equation, substitution for Ci
        Farquhar, G.D., and S. von Caemmerer, 1982.  Modelling of photosynthetic
                response to environmental conditions.  In Encyclopedia of Plant
                Physiology, New Series, Vol. 12B, Physiological Plant Ecology II,
                O.L. Lange, P.S. Nobel, C.B. Osmond, and H. Ziegler, eds, Springer-
                Verlag, Berlin, Germany, pp 549-587.

        from A = g(Ca-Ci) into the equations from Farquhar and von Caemmerer:

               Vcmax (Ci - gammaStar)
        Av =  -------------------   -   Rd
              Ci + Kc (1 + O2/Ko)

        Use Aj equation from dePury and Farquhar, 1997

                 J (Ci - gammaStaR)
        Aj  =  -------------------  -   Rd
              4(Ci +  2*gammaStar)

   	and As as sink-limited photosynthesis from accumulation of starch at the chloroplast
                (Bonan 2008, page 246)

        As = Vxmax / 2

        */

// quadratic solution for Av
        aa = -1.0/g;
        bb = Ca + (Vcmax - Rd)/g + Kc*(1.0 + O2/Ko);
        cc = Vcmax*(gammaStar - Ca) + Rd*(Ca + Kc*(1.0 + O2/Ko));

        if ((det = bb*bb - 4.0*aa*cc) < 0.0)
        {
                Av = 0.0;
        }
        else
        {
                Av = (-bb + sqrt(det)) / (2.0*aa);
        }

// quadratic solution for Aj
        aa = -4.0/g;
        bb = 4.0*Ca + 8.0*gammaStar + J/g - 4.0*Rd/g;
        cc = J*(gammaStar - Ca) + Rd*(4.0*Ca + 8.0*gammaStar);

        if ((det = bb*bb - 4.0*aa*cc) < 0.0)
        {
                Aj = 0.0;
        }
        else
        {
                Aj = (-bb + sqrt(det)) / (2.0*aa);
        }
//cout << "Aj = " << Aj << endl;
// sink-limited As
        As = Vcmax / 2.0;

//determine Rubisco activity or electron transport limited
        if (Av < Aj)
        {
                alpha1 = Vcmax;
                alpha2 = Kc*(1.0+O2/Ko);
                g_mult = g;
        }
        else
        {
                alpha1 = J;
                alpha2 = 2.0*gammaStar;
                g_mult = 4.0*g;
        }

//implementation of quadratic solution of Ci
//from Katul et al, 2003, Plant, Cell & Environment
//DSM June 7, 2005 - Added Rd
        kappa1 = -g_mult;
        kappa2 = g_mult * (Ca - alpha2) - alpha1;
        kappa3 = g_mult*alpha2*Ca + alpha1*gammaStar;

        Ci = (-kappa2 - sqrt(kappa2*kappa2 - 4.0*kappa1*kappa3))/(2.0*kappa1);

        if (Ci > (Ca+0.000001))
        {
                Ci = (-kappa2 + sqrt(kappa2*kappa2 - 4.0*kappa1*kappa3))/(2.0*kappa1);
        }
        if (Ci < gammaStar)
        {
                Ci = gammaStar;
        }

    //Calculate conductance using most limiting of Rubisco-limited or
    //electron-transport limited assimilation, and Ci from Katul
    
    if (Av > Aj)
    {
        A = Aj;
    }
    else
    {
        A = Av;
    }
    
    if (A > As)
      {
             A = As;
    }
    
    A_kg = A / 4.6296; //convert from umol m-2 s-1 to kg ha-1
    
    //Positive photosynthesis
    if (A > 0.0)
    {
        g = (A)/(Ca-Ci+0.01);
        out.g = g * (R * tk) / 1.0e6;
    }
    //At compensation point or lower, hold stomata open
    else
    {
        if (getChloroplastStarch() > (-A_kg))
        {
            putChloroplastStarch(getChloroplastStarch()-(-A_kg));
        }
        else
        {
            Rd = Rd * getChloroplastStarch()/(-A_kg);
            A = A -getChloroplastStarch()*4.6296;
            putChloroplastStarch(0.0);
        }
        g = (-A)/(2.0*gammaStar);
        out.g = g * (R * tk) / 1.0e6;
    }
    //Cuticular conductance
        if (g < 0.0001/1.6/42.0)
        {
                g = 0.0001/1.6/42.0;
                out.g = g;
        }
        D = in.D;
        E = 1.6*g*D;

        out.E = E;
        out.Ca = Ca; 			//(Pa) atmospheric [CO2]
        out.Ci = Ci*1.0e6/in.pa;	//(Pa) intercellular [CO2]
        out.gammaStar = gammaStar;	//(Pa) CO2 compensation point, no Rd
        out.O2 = O2;			//(Pa) atmospheric [O2]
        out.Kc = Kc;			//(Pa) MM constant carboxylation
        out.Ko = Ko;			//(Pa) MM constant oxygenation
        out.act = act;			//(umol/kg/s) Rubisco activity
        out.Vcmax25 = Vcmax25;		//(umol/m2/s) max rate carboxylation at 25C
        out.Vcmax = Vcmax;		//(umol/m2/s) max rate carboxylation
        out.Jmax25 = Jmax25;		//(umol/m2/s) max rate electron transport at 25C
        out.Jmax = Jmax;		//(umol/m2/s) max rate electron transport
        out.J = J;			//(umol/m2/s) rate of RuBP regeneration
        out.Av = Av;			//(umol/m2/s) carboxylation limited assimilation
        out.Aj = Aj;			//(umol/m2/s) RuBP regen limited assimilation
        out.A = A;			//(umol/m2/s) final assimilation rate
        out.betaA = betaA;
        out.phi2 = phi2;
        out.s = s;
        out.alpha = alpha;
        return out.A;
} //end photosynthesis function


//
//computeRootArea()
//calculate absorbing root area based on the fine root carbon
//Limited data suggests a SRL average in the first two root classes of about
//    30 m gC-1 in conifer and about twice as much in deciduous
//    e.g., Withington et al 2006 Ecological Monographs
//		Meinen et al 2009 Oecologia
//
double BiogeochemicalCycles::computeRootArea(trees_params treesParams)
{
        double totalRootArea;

        totalRootArea = 0.0;

        for (int j = 0; j < nRoots; j++)
        {
		totalRootArea += computeRootArea(treesParams, j);
        }
	return(totalRootArea);
}
double BiogeochemicalCycles::computeRootArea(trees_params treesParams,
					     int j)
{
        double rootCarbon, totalRootArea, SRL, SRA, rootDiam, refRootDiam;

	assert(j >= 0);
	assert(j < nRoots);

        totalRootArea = 0.0;
	refRootDiam = 0.00025; //m
	//rootDiam = treesParams.maxRootDiam * pow(0.5, 9);
	rootDiam = 0.000125; //m
//specific root length (m kgC-1 root)
//SRL1 is defined at a root diameter of 0.00025 m
        SRL = treesParams.SRL1*pow(refRootDiam, 2.0)/pow(rootDiam,2.0);
        for (int k = 0; k < nRootOrders; k++)
        {
//root segments are approximated as cylinders
//specific root area = m2 root m-2 ground area kgC-1 root
                SRA = SRL * rootDiam * M_PI;
		rootCarbon = fineRootBiomassCarbon[j][k] + coarseRootBiomassCarbon[j][k];
//convert kgC ha-1 to kgC m-2
//m2 root area m-2 ground area = 
//      kgC root m-2 ground area * m-2 ground area * m2 root m-2 ground area kgC-1 root
		rootArea[j][k] = 0.0001 * rootCarbon * SRA;
                totalRootArea += rootArea[j][k];
                rootDiam *= treesParams.rootDiamMultiplier;
                SRL /= pow(treesParams.rootDiamMultiplier, 2.0);
        }
	return(totalRootArea);
}

//
//computeFineRootArea()
//calculate absorbing root area based on the fine root carbon only
//Limited data suggests a SRL average in the first two root classes of about
//    30 m gC-1 in conifer and about twice as much in deciduous
//    e.g., Withington et al 2006 Ecological Monographs
//		Meinen et al 2009 Oecologia
//
double BiogeochemicalCycles::computeFineRootArea(trees_params treesParams)
{
        double totalRootArea;

        totalRootArea = 0.0;
        for (int j = 0; j < nRoots; j++)
        {
		totalRootArea += computeFineRootArea(treesParams, j);
        }
	return(totalRootArea);
}
double BiogeochemicalCycles::computeFineRootArea(trees_params treesParams,
					         int j)
{
        double rootCarbon, totalRootArea, SRL, SRA, rootDiam, refRootDiam;

	assert(j >= 0);
	assert(j < nRoots);

        totalRootArea = 0.0;
	refRootDiam = 0.00025; //m
	//rootDiam = treesParams.maxRootDiam * pow(0.5, 9);
	rootDiam = 0.000125; //m
//specific root length (m kgC-1 root)
//SRL1 is defined at a root diameter of 0.00025 m
        SRL = treesParams.SRL1*pow(refRootDiam, 2.0)/pow(rootDiam,2.0);
        for (int k = 0; k < nFineRootOrders; k++)
        {
//specific root length (m kgC-1 root)
//root segments are approximated as cylinders
                SRA = SRL * rootDiam * M_PI;
		rootCarbon = fineRootBiomassCarbon[j][k] + coarseRootBiomassCarbon[j][k];
//convert kgC ha-1 to kgC m-2
//m2 root area m-2 ground area = 
//      kgC root m-2 ground area * m-2 ground area * m2 root m-2 ground area kgC-1 root
		rootArea[j][k] = 0.0001 * rootCarbon * SRA; //m2 root m2 ground area
                totalRootArea += rootArea[j][k];
                rootDiam *= treesParams.rootDiamMultiplier;
                SRL /= pow(treesParams.rootDiamMultiplier, 2.0);
                //SRL /= 4.0;
        }
	return(totalRootArea);
}

//
//computeLateralRootWeibulls(), updateLateralRootWeibulls()
//  This code adjusts the lateral root saturated vulnerability curve parameters.
//  Adjustment of b is downward for increases in fine root area and increased
//    fine root CN ratio. The adjustment of c is upward for the same conditions.
//
void BiogeochemicalCycles::computeLateralRootWeibulls(trees_params& treesParams)
{
        double bval, cval, bsum, csum;
        double totalRootArea, totalRootCarbon, rootCarbonRatio;

        totalRootArea = computeRootArea(treesParams);
        for (int j = 0; j < nRoots; j++)
        {
        	bsum = csum = 0.0;
                bval = lat_Root_b_value_init[0];
                cval = lat_Root_c_value_init[0];
		totalRootCarbon = getFineRootCarbon(j);
		rootCarbonRatio = fineRootBiomassCarbon[j][2]/totalRootCarbon;
                bsum += bval*rootCarbonRatio;
                csum += cval*rootCarbonRatio;
                for (int k = 1; k >= 0; k--)
                {
                        bval *= 0.95;
                        cval /= 0.95;
			rootCarbonRatio = fineRootBiomassCarbon[j][k]/totalRootCarbon;
                        bsum += bval*rootCarbonRatio;
                        csum += cval*rootCarbonRatio;
                }
                bval = lat_Root_b_value_init[0];
                cval = lat_Root_c_value_init[0];
                for (int k = 3; k <= 4; k++)
                {
                        bval /= 0.95;
                        cval *= 0.95;
			rootCarbonRatio = fineRootBiomassCarbon[j][k]/totalRootCarbon;
                        bsum += bval*rootCarbonRatio;
                        csum += cval*rootCarbonRatio;
                }
        }
        treesParams.lat_Root_b_value = bsum;
        treesParams.lat_Root_c_value = csum;
}
void BiogeochemicalCycles::computeLateralRootWeibulls(int j, 
							trees_params treesParams,
							double& bsum,
							double& csum)
{
        double bval, cval;
        double totalRootArea, CN, CNratio, SLAscalar, rootScalar;
	double totalRootCarbon, rootCarbonRatio;
	double bfactor = 0.95;
	double cfactor = 1.0/0.95;
	SLAscalar = 22.0/treesParams.SLA;
	if (SLAscalar > 3.0)
	{
		SLAscalar = 3.0;
	}
	else if (SLAscalar < 1.0)
	{
		SLAscalar = 1.0;
	}
	rootScalar = 1.5625;
        bsum = csum = 0.0;
	totalRootCarbon = getFineRootCarbon(j);
	CN = fineRootBiomassCarbon[j][2]/fineRootBiomassNitrogen[j][2];
	CN = max(CN, 20.0);
	CNratio = 0.5 * (1.0 + 20.0*rootScalar*SLAscalar/CN);
	CNratio = max(0.5, CNratio);
        bval = lat_Root_b_value_init[0]*CNratio;
        cval = lat_Root_c_value_init[0]/CNratio;
	rootCarbonRatio = fineRootBiomassCarbon[j][2]/totalRootCarbon;
        bsum += bval*rootCarbonRatio;
        csum += cval*rootCarbonRatio;
        for (int k = 1; k >= 0; k--)
        {
		rootScalar /= 1.25;
		CN = fineRootBiomassCarbon[j][k]/fineRootBiomassNitrogen[j][k];
		CN = max(CN, 20.0);
		CNratio = 0.5 * (1.0 + 20.0*rootScalar*SLAscalar/CN);
		CNratio = max(0.5, CNratio);
        	bval = lat_Root_b_value_init[0]*(bfactor*CNratio);
                cval = lat_Root_c_value_init[0]*(cfactor/CNratio);
		rootCarbonRatio = fineRootBiomassCarbon[j][k]/totalRootCarbon;
                bsum += bval*rootCarbonRatio;
                csum += cval*rootCarbonRatio;
		bfactor *= 0.95;
		cfactor /= 0.95;
        }
	rootScalar = 1.5625;
	bfactor = 1.0/0.95;
	cfactor = 0.95;
        bval = lat_Root_b_value_init[0]*bfactor;
        cval = lat_Root_c_value_init[0]*cfactor;
        for (int k = 3; k <= 4; k++)
        {
		rootScalar *= 1.25;
		CN = fineRootBiomassCarbon[j][k]/fineRootBiomassNitrogen[j][k];
		CN = max(CN, 20.0);
		CNratio = 20.0*rootScalar*SLAscalar/CN;
        	bval = lat_Root_b_value_init[0]*(bfactor*CNratio);
                cval = lat_Root_c_value_init[0]*(cfactor/CNratio);
		rootCarbonRatio = fineRootBiomassCarbon[j][k]/totalRootCarbon;
                bsum += bval*rootCarbonRatio;
                csum += cval*rootCarbonRatio;
		bfactor /= 0.95;
		cfactor *= 0.95;
        }
}
void BiogeochemicalCycles::updateLateralRootWeibulls(double bsat[][MD], 
							double ccsat[][MD],
							trees_params treesParams)
{
	double bsum, csum;
        for (int j = treesParams.smodules+2; j < nRoots+treesParams.smodules+2; j++)
        {
        	bsum = csum = 0.0;
		computeLateralRootWeibulls(j-treesParams.smodules-2, treesParams, bsum, csum);
		bsat[j][1] = bsum;
		ccsat[j][1] = csum;
        }
}

//
//getLeafBiomassCarbon()
//
double BiogeochemicalCycles::getLeafBiomassCarbon()
{
	return(leafBiomassCarbon[0]);
}

//
//getRootCarbon()
//This code sums up all the carbon within the root system
//
double BiogeochemicalCycles::getRootCarbon()
{
	double rootCarbon = 0.0;

	for (int j = 0; j < nRoots; j++)
	{
		for (int k = 0; k < nRootOrders; k++)
		{
			rootCarbon += getRootCarbon(j, k);
		}
	}
	return (rootCarbon);
}
double BiogeochemicalCycles::getRootCarbon(int j, 
					   int k)
{
	double rootCarbon;

	assert(j >= 0);
	assert(j < nRoots);
	assert(k >= 0);
	assert(k < nRootOrders);

	rootCarbon = fineRootBiomassCarbon[j][k] + coarseRootBiomassCarbon[j][k];
	return(rootCarbon);
}

//
//getFineRootCarbon()
//This code sums up all the fine root carbon within the root system
// --Note: Assumption is that the first 5 root orders are fine (0.25 mm to 4.0 mm diameter)
//
double BiogeochemicalCycles::getFineRootCarbon()
{
	double rootCarbon = 0.0;

	for (int j = 0; j < nRoots; j++)
	{
		rootCarbon += getFineRootCarbon(j);
	}
	return (rootCarbon);
}
double BiogeochemicalCycles::getFineRootCarbon(int j)
{
	double rootCarbon = 0.0;

	for (int k = 0; k < nFineRootOrders; k++)
	{
		rootCarbon += getRootCarbon(j, k);
	}
	return (rootCarbon);
}

//
//getLiveStemCarbon()
//
double BiogeochemicalCycles::getLiveStemCarbon()
{
	double stemCarbon;
	stemCarbon = liveStemCarbon[0];
	return(stemCarbon);
}

//
//getDeadStemCarbon()
//
double BiogeochemicalCycles::getDeadStemCarbon()
{
	double stemCarbon;
	stemCarbon = deadStemCarbon[0];
	return(stemCarbon);
}

//
//getLeafNSC()
//
double BiogeochemicalCycles::getLeafNSC()
{
	double nsc;
	nsc = leafNSC[0];
	return(nsc);
}

//
//getChloroplastStarch()
//
double BiogeochemicalCycles::getChloroplastStarch()
{
	double nsc;
	nsc = chloroplastStarch[0];
	return(nsc);
}

//
//getChloroplastSugar()
//
double BiogeochemicalCycles::getChloroplastSugar()
{
	double nsc;
	nsc = chloroplastSugar[0];
	return(nsc);
}

//
//getStemNSC()
//
double BiogeochemicalCycles::getStemNSC()
{
	double nsc;
	nsc = stemNSC[0];
	return(nsc);
}

//
//getLeafN()
//
double BiogeochemicalCycles::getLeafBiomassN()
{
        return(leafBiomassNitrogen[0]);
}

//
//getRootNSC()
//
double BiogeochemicalCycles::getRootNSC()
{
	double nsc = 0.0;

	for (int j = 0; j < nRoots; j++)
	{
		for (int k = 0; k < nRootOrders; k++)
		{
			nsc += getRootNSC(j, k);
		}
	}
	return(nsc);
}

//
//getRootNSC()
//
double BiogeochemicalCycles::getRootNSC(int j,
					int k)
{
	double nsc;

	assert(j >= 0);
	assert(j < nRoots);
	assert(k >= 0);
	assert(k < nRootOrders);

	nsc = rootNSC[j][k];
	return(nsc);
}

//
//putLeafNSC()
//
void BiogeochemicalCycles::putLeafNSC(double nsc)
{
	leafNSC[0] = nsc;
}

//
//putChloroplastStarch()
//
void BiogeochemicalCycles::putChloroplastStarch(double nsc)
{
	chloroplastStarch[0] = nsc;
}

//
//putChloroplastSugar()
//
void BiogeochemicalCycles::putChloroplastSugar(double nsc)
{
	chloroplastSugar[0] = nsc;
}

//
//putStemNSC()
//
void BiogeochemicalCycles::putStemNSC(double nsc)
{
	stemNSC[0] = nsc;
}

//
//updateRootNSC()
//This functions serves primiarly as a means of deducting NSC from roots
//   as a cost of root respiration
//Use the first function if you don't have root-specific NSC updates, as this
//   function will distribute the NSC consumption in proportion to the live
//   root carbon.
//The better approach is to have root-specific NSC updates, in which case
//   you would use the second function.
//
void BiogeochemicalCycles::updateRootNSC(double delta_nsc)
{
	double rootC, rootCproportion;

	rootC = getRootCarbon();
	if (rootC > 0.00001)
	{
		for (int j = 0; j < nRoots; j++)
		{
			for (int k = 0; k < nRootOrders; k++)
			{
				rootCproportion = getRootCarbon(j,k)/rootC;
				rootNSC[j][k] += delta_nsc*rootCproportion;
				if (rootNSC[j][k] < 0.000001)
				{
					rootNSC[j][k] = 0.000001;
				}
			}
		}
	}
}
void BiogeochemicalCycles::updateRootNSC(double delta_nsc,
					 int j,
					 int k)
{

	assert(j >= 0);
	assert(j < nRoots);
	assert(k >= 0);
	assert(k < nRootOrders);

	rootNSC[j][k] += delta_nsc;

}

//
//getHumus()
//
double BiogeochemicalCycles::getHumus()
{
	double humus = 0.0;

	for (int j = 0; j < nRoots; j++)
	{
		humus += getHumus(j);
	}
	return(humus);
}
double BiogeochemicalCycles::getHumus(int j)
{
	double humus;

	assert(j >= 0);
	assert(j < nRoots);

	humus = humusCarbon[j];
	return humus;
}

//
//getRhizosphereCl()
//
double BiogeochemicalCycles::getRhizosphereCl()
{
	double Cl = 0.0;
	for (int j = 0; j < nRoots; j++)
	{
		for (int k = 0; k < nRootOrders; k++)
		{
			Cl += getRhizosphereCl(j, k);
		}
	}
	return(Cl);
}
double BiogeochemicalCycles::getRhizosphereCl(int j)
{
	double Cl = 0.0;

	assert(j >= 0);
	assert(j < nRoots);

	for (int k = 0; k < nRootOrders; k++)
	{
		Cl += getRhizosphereCl(j, k);
	}
	return(Cl);
}
double BiogeochemicalCycles::getRhizosphereCl(int j,
 					      int k)
{
	double Cl;

	assert(j >= 0);
	assert(j < nRoots);
	assert(k >= 0);
	assert(k < nRootOrders);

	Cl = rhizosphereCl[j][k];
	return(Cl);
}

//
//getRhizosphereNl()
//
double BiogeochemicalCycles::getRhizosphereNl()
{
	double Nl = 0.0;
	for (int j = 0; j < nRoots; j++)
	{
		for (int k = 0; k < nRootOrders; k++)
		{
			Nl += getRhizosphereNl(j, k);
		}
	}
	return(Nl);
}
double BiogeochemicalCycles::getRhizosphereNl(int j)
{
	double Nl = 0.0;

	assert(j >= 0);
	assert(j < nRoots);

	for (int k = 0; k < nRootOrders; k++)
	{
		Nl += getRhizosphereNl(j, k);
	}
	return(Nl);
}
double BiogeochemicalCycles::getRhizosphereNl(int j,
					      int k)
{
	double Nl;

	assert(j >= 0);
	assert(j < nRoots);
	assert(k >= 0);
	assert(k < nRootOrders);

	Nl = rhizosphereNl[j][k];
	return(Nl);
}

//
//getRhizosphereLiveMicrobialCarbon()
//
double BiogeochemicalCycles::getRhizosphereLiveMicrobialCarbon()
{
        double Cm = 0.0;
        for (int j = 0; j < nRoots; j++)
        {
                for (int k = 0; k < nRootOrders; k++)
                {
                        Cm += getRhizosphereLiveMicrobialCarbon(j, k);
                }
        }
        return(Cm);
}
double BiogeochemicalCycles::getRhizosphereLiveMicrobialCarbon(int j)
{
	double Cm = 0.0;

	assert(j >= 0);
	assert(j < nRoots);

	for (int k = 0; k < nRootOrders; k++)
	{
		Cm += getRhizosphereLiveMicrobialCarbon(j, k);
	}
	return(Cm);
}
double BiogeochemicalCycles::getRhizosphereLiveMicrobialCarbon(int j,
                                            		       int k)
{
        double Cm;

	assert(j >= 0);
	assert(j < nRoots);
	assert(k >= 0);
	assert(k < nRootOrders);

        Cm = rhizosphereLiveMicrobialCarbon[j][k];
        return(Cm);
}

//
//getRhizosphereMicrobialNitrogen()
//
double BiogeochemicalCycles::getRhizosphereMicrobialNitrogen()
{
        double Nm = 0.0;
        for (int j = 0; j < nRoots; j++)
        {
                for (int k = 0; k < nRootOrders; k++)
                {
                        Nm += getRhizosphereMicrobialNitrogen(j, k);
                }
        }
        return(Nm);
}
double BiogeochemicalCycles::getRhizosphereMicrobialNitrogen(int j)
{
	double Nm = 0.0;

	assert(j >= 0);
	assert(j < nRoots);

	for (int k = 0; k < nRootOrders; k++)
	{
		Nm += getRhizosphereMicrobialNitrogen(j, k);
	}
	return(Nm);
}
double BiogeochemicalCycles::getRhizosphereMicrobialNitrogen(int j,
                                            		     int k)
{
        double Nm;

	assert(j >= 0);
	assert(j < nRoots);
	assert(k >= 0);
	assert(k < nRootOrders);

        Nm = rhizosphereMicrobialNitrogen[j][k];
        return(Nm);
}

//
//getRhizosphereDeadMicrobialCarbon()
//
double BiogeochemicalCycles::getRhizosphereDeadMicrobialCarbon()
{
        double Cm = 0.0;
        for (int j = 0; j < nRoots; j++)
        {                                   
                for (int k = 0; k < nRootOrders; k++)
                {
                        Cm += getRhizosphereDeadMicrobialCarbon(j, k);
                }
        }
        return(Cm);
}
double BiogeochemicalCycles::getRhizosphereDeadMicrobialCarbon(int j,
                                            		       int k)
{
        double Cm;

	assert(j >= 0);
	assert(j < nRoots);
	assert(k >= 0);
	assert(k < nRootOrders);

        Cm = rhizosphereDeadMicrobialCarbon[j][k];
        return(Cm);
}

//
//getRhizosphereMineralNitrogen()
//
double BiogeochemicalCycles::getRhizosphereMineralNitrogen()
{
        double Nm = 0.0;
        for (int j = 0; j < nRoots; j++)
        {
                for (int k = 0; k < nRootOrders; k++)
                {
                        Nm += getRhizosphereMineralNitrogen(j, k);
                }
        }
        return(Nm);
}
double BiogeochemicalCycles::getRhizosphereMineralNitrogen(int j)
{
        double Nm = 0.0;
        for (int k = 0; k < nRootOrders; k++)
        {
                Nm += getRhizosphereMineralNitrogen(j, k);
        }
        return(Nm);
}
double BiogeochemicalCycles::getRhizosphereMineralNitrogen(int j,
                                            		   int k)
{
        double Nm;
        Nm = rhizosphereMineralNitrogen[j][k];
        return(Nm);
}

double BiogeochemicalCycles::getNitrogenLeaching(int j)
{
	double Nl;
	Nl = nitrogenLeaching[j];
	return(Nl);
}

//
//getRhizosphereAmmoniumNitrogen()
//
double BiogeochemicalCycles::getRhizosphereAmmoniumNitrogen()
{
        double Nm = 0.0;
        for (int j = 0; j < nRoots; j++)
        {
                for (int k = 0; k < nRootOrders; k++)
                {
                        Nm += getRhizosphereAmmoniumNitrogen(j, k);
                }
        }
        return(Nm);
}
double BiogeochemicalCycles::getRhizosphereAmmoniumNitrogen(int j)
{
        double Nm = 0.0;
        for (int k = 0; k < nRootOrders; k++)
        {
                Nm += getRhizosphereAmmoniumNitrogen(j, k);
        }
        return(Nm);
}
double BiogeochemicalCycles::getRhizosphereAmmoniumNitrogen(int j,
                                            		   int k)
{
        double Nm;
        Nm = rhizosphereAmmoniumNitrogen[j][k];
        return(Nm);
}

//
//getRhizosphereNitrateNitrogen()
//
double BiogeochemicalCycles::getRhizosphereNitrateNitrogen()
{
        double Nm = 0.0;
        for (int j = 0; j < nRoots; j++)
        {
                for (int k = 0; k < nRootOrders; k++)
                {
                        Nm += getRhizosphereNitrateNitrogen(j, k);
                }
        }
        return(Nm);
}
double BiogeochemicalCycles::getRhizosphereNitrateNitrogen(int j)
{
        double Nm = 0.0;
        for (int k = 0; k < nRootOrders; k++)
        {
                Nm += getRhizosphereNitrateNitrogen(j, k);
        }
        return(Nm);
}
double BiogeochemicalCycles::getRhizosphereNitrateNitrogen(int j,
                                            		   int k)
{
        double Nm;
        Nm = rhizosphereNitrateNitrogen[j][k];
        return(Nm);
}

//
//computeNSCfluxes(kratio, kratio_vector)
//Move NSC between leaf-stem and stem-root pools based on concentration gradients
//  and relative hydraulic conductance (kratio)
//
//Assumes phloem loading has already taken place and energy costs consumed
//  -- these costs are taken where net photosynthesis is added to the PSN state variable in 
//     simulation_function()
//
//Other assumptions: 0.08 diffusion term assumes that it would take 12 simulation time
//                   steps or 6 hours to completely move the NSC from leaf to root
//
void BiogeochemicalCycles::computeNSCfluxes(double kratio,
				       	    double* kratio_vector)
{
	double C, flux, leafCfract, stemCfract, rootCfract; 
	double rootCarbon, singleRootC, kratioRoot;
	double cfract;
//flux between stem and roots
	rootCarbon = getRootCarbon();
	for (int j = 0; j < nRoots; j++)
	{
		kratioRoot = kratio_vector[j];
//nFineRootOrders
		cfract = 0.1;
		for (int k = 0; k < nRootOrders; k++)
		{
			singleRootC = fineRootBiomassCarbon[j][k]+coarseRootBiomassCarbon[j][k];
			rootCfract = rootNSC[j][k]/(singleRootC+0.0000001);
			stemCfract = stemNSC[0]/liveStemCarbon[0];
			if (stemCfract > cfract*rootCfract)
			{
				C = 0.20*stemNSC[0]*(singleRootC/rootCarbon)*kratioRoot;
				//C = 0.020*stemNSC[0]*kratioRoot;
			}
			else if (stemCfract < cfract*rootCfract)
			{
				C = 0.08*rootNSC[j][k]*(singleRootC/rootCarbon)*kratioRoot;
			}
			else
			{
				C = 0.0;
			}
			flux = C*(stemCfract - cfract*rootCfract);
			rootNSC[j][k] += flux;
			stemNSC[0] -= flux;
			cfract += 0.1;
		}
	}
//flux between canopy and stem
	leafCfract = leafNSC[0]/leafBiomassCarbon[0];
	stemCfract = stemNSC[0]/liveStemCarbon[0];
	if (leafCfract > 1.0*stemCfract)
	{
		C = 0.08*leafNSC[0]*kratio;
		//C *= leafCfract/0.1;
	}
	else if (leafCfract < 1.0*stemCfract)
	{
		C = 0.20*stemNSC[0]*kratio;
	}
	else
	{
		C = 0.0;
	}
	flux = C*(leafCfract-1.0*stemCfract);
	leafNSC[0] -= flux;
	stemNSC[0] += flux;
}


void BiogeochemicalCycles::computeNSCfluxes(double kratio,
                                            double kratioRoot)
{
        double C, flux, leafCfract, stemCfract, rootCfract, rootCarbon, singleRootC;
	double cfract;
//flux between stem and roots
        rootCarbon = getRootCarbon();
        for (int j = 0; j < nRoots; j++)
        {
		cfract = 0.1;
                for (int k = 0; k < nRootOrders; k++)
                {
                        singleRootC = fineRootBiomassCarbon[j][k] + coarseRootBiomassCarbon[j][k];
                        rootCfract = rootNSC[j][k]/singleRootC;
                        stemCfract = stemNSC[0]/liveStemCarbon[0];
                        if (stemCfract > cfract*rootCfract)
                        {
                                C = 0.20*stemNSC[0]*(singleRootC/rootCarbon)*kratioRoot;
                                //C = 0.020*stemNSC[0]*kratioRoot;
                        }
                        else if (stemCfract < cfract*rootCfract)
                        {
                                C = 0.08*rootNSC[j][k]*(singleRootC/rootCarbon)*kratioRoot;
                        }
                        else
                        {
                                C = 0.0;
                        }
                        flux = C*(stemCfract - cfract*rootCfract);
                        rootNSC[j][k] += flux;
                        stemNSC[0] -= flux;
			cfract += 0.1;
                }
        }
//flux between canopy and stem
        leafCfract = leafNSC[0]/leafBiomassCarbon[0];
        stemCfract = stemNSC[0]/liveStemCarbon[0];
        if (leafCfract > 2.5*stemCfract)
        {
                C = 0.08*leafNSC[0]*kratio;
        }
        else if (leafCfract < 2.5*stemCfract)
        {
                C = 0.20*stemNSC[0]*kratio;
        }
        else
        {
                C = 0.0;
        }
        flux = C*(leafCfract-2.5*stemCfract);
        leafNSC[0] -= flux;
        stemNSC[0] += flux;
}


//
//computeNfluxes(kratio, kratio_vector)
//Move N between leaf-root pools based on concentration gradients
//    and relative hydraulic conductance
//
void BiogeochemicalCycles::computeNfluxes(double kratio,
					  double* kratio_vector)
{
	double C, cfract, flux, leafNfract, rootNfract, rootNitrogen, leafNstore, singleRootN, kratioRoot;

	leafNstore = leafStoredNitrogen[0];
	rootNitrogen = getRootN();
	for (int j = 0; j < nRoots; j++)
	{
		cfract = 0.1;
		kratioRoot = kratio_vector[j];
		for (int k = 0; k < nRootOrders; k++)
		{
			singleRootN = fineRootBiomassNitrogen[j][k] + coarseRootBiomassNitrogen[j][k];
			leafNfract = leafStoredNitrogen[0]/leafBiomassNitrogen[0];
			rootNfract = rootMineralNitrogen[j][k]/(singleRootN+0.0000001);
			if (leafNfract > cfract*rootNfract)
			{
				C = 0.08*leafNstore*(singleRootN/rootNitrogen)*kratioRoot;
				//C = 0.016*leafNstore*kratioRoot;
			}
			else if (leafNfract < cfract*rootNfract)
			{
				C = 0.08*rootMineralNitrogen[j][k]*(singleRootN/rootNitrogen)*kratioRoot;
			}
			else
			{
				C = 0.0;
			}
			flux = C*(leafNfract - cfract*rootNfract);
			rootMineralNitrogen[j][k] += flux;
			leafStoredNitrogen[0] -= flux;
		}
	}
}

void BiogeochemicalCycles::computeNfluxes(double kratio,
                                          double kratioRoot)
{
        double C, flux, leafNfract, rootNfract, rootNitrogen, leafNstore, singleRootN;

        leafNfract = leafStoredNitrogen[0]/leafBiomassNitrogen[0];
        leafNstore = leafStoredNitrogen[0];
        rootNitrogen = getRootN();
        for (int j = 0; j < nRoots; j++)
        {
                for (int k = 0; k < nRootOrders; k++)
                {
                        singleRootN = fineRootBiomassNitrogen[j][k] + coarseRootBiomassNitrogen[j][k];
                        rootNfract = rootMineralNitrogen[j][k]/singleRootN;
                        if (leafNfract > 0.5*rootNfract)
                        {
                                C = 0.08*leafNstore*(singleRootN/rootNitrogen)*(1.0/(1.0/kratio + 1.0/kratioRoot));
                                //C = 0.016*leafNstore*(1.0/(1.0/kratio + 1.0/kratioRoot));
                        }
                        else if (leafNfract < 0.5*rootNfract)
                        {
                                C = 0.08*rootMineralNitrogen[j][k]*(singleRootN/rootNitrogen)*(1.0/(1.0/kratio + 1.0/kratioRoot));
                        }
                        else
                        {
                                C = 0.0;
                        }
                        flux = C*(leafNfract - 1.0*rootNfract);
                        rootMineralNitrogen[j][k] += flux;
                        leafStoredNitrogen[0] -= flux;
                }
        }
}


//
//computeLeafAllocation()
//Update growth and leaf area and root-to-leaf area ratio if not in dormancy (root T < 5 C)
//What to do about LAI - need to modify lai and Al; lai_at_sat_kl is used only once per simulation
//Currently assumes allocation to new leaf growth is up to 45% of growth respiration, 
//    and leaf longevity is an input
//Assuming leaf biomass can be approximated by 86% (carbon in cellulose) of NSC use
//This computes potential lai at saturated kl, with actual determined by phenology
//NEED - LIVE LAI for canopy conductance since brown canopy LAI absorbs energy, but does not transpire
//
void BiogeochemicalCycles::computeLeafAllocation(trees_params& treesParams,
						 double newc[][MD],
						 double N_avail_rate_plant,
						 double kratio,
						 double nsc,
						 double psn,
						 double nscRatio,
						 double& rgrowth,
						 double& leafCfraction,
						 double lai,
						 double& stressedLeafLifeSpan,
						 double SLA)
{
	double deltaLAI, unstressedLeafLifeSpan;

	deltaLAI = 0.0;
        unstressedLeafLifeSpan = treesParams.leafLifeSpan*365.25*48;
//start with a base fraction of carbon allocation to leaf
        leafCfraction = 0.40;
//if available N is low then reduce allocation to leaf
	if (N_avail_rate_plant < 0.5)
	{
		N_avail_rate_plant *= 2.0;
	}
	else
	{
		N_avail_rate_plant = 1.0;
	}
        leafCfraction = leafCfraction*N_avail_rate_plant*kratio;
        //leafCfraction = leafCfraction*N_avail_rate_plant;

        if (leafCfraction < 0.0001)
        {
                leafCfraction = 0.0001;
        }
        else if (leafCfraction > 0.40)
        {
                leafCfraction = 0.40;
        }

//cout << "leafCfraction = " << leafCfraction << endl;

	nsc = getLeafNSC()+getStemNSC()+getRootNSC();

//prevent growth allocation of bringing NSC to zero
//increase or decrease growth by a nonlinear function of excess or deficit NSC
	if (rgrowth > 0.99*nsc)
	{
		rgrowth = 0.99*nsc;
	}

//lifespan in units of 30 minutes
//compute unstressed change in LAI
        deltaLAI = leafCfraction*rgrowth*treesParams.SLA/10000.0 - treesParams.lai_at_sat_kl/unstressedLeafLifeSpan;
        treesParams.lai_at_sat_kl += deltaLAI;

	if (kratio < 0.001)
	{
		kratio = 0.001;
	}
//set the leaf stress level if it has K < 50% Ksat and stress is higher in leaf than in shallow root
        if (newc[1][1] > newc[3][1])
        {
                stressedLeafLifeSpan = kratio*unstressedLeafLifeSpan;
        }
        else
        {
                stressedLeafLifeSpan = unstressedLeafLifeSpan;
        }
//compute stressed change in LAI
        deltaLAI = leafCfraction*rgrowth*treesParams.SLA/10000.0 - treesParams.live_lai/stressedLeafLifeSpan;
        treesParams.live_lai += deltaLAI;

        if (treesParams.lai < 0.01)
        {
                treesParams.lai = 0.01;
        }
        if (treesParams.lai_at_sat_kl < 0.01)
        {
                treesParams.lai_at_sat_kl = 0.01;
        }
        if (treesParams.live_lai < 0.01)
        {
                treesParams.live_lai = 0.01;
        }
}

//
//updateLeafCarbonNitrogenPools()
//
void BiogeochemicalCycles::updateLeafCarbonNitrogenPools(trees_params& treesParams,
						 	 double delta_lai,
						 	 double RL,
							 double N_neg_fract,
							 double& N_neg_demand,
							 double& N_pos_demand)
{
	double CN, minLeafC, baseLAI, nsc, delta_nsc, GRcost;
	double total_N_demand, residual_N_demand, NfromStorage;
	double N_avail_rate_plant = plantNstatus[0];
//when delta_lai is positive, then increase leaf biomass
//when delta_lai is negative, then retranslate 50% of leaf protein N and reduce biomass

//set leaf carbon minimum for the leaf bud to equivalent of 0.1 m2 m-2 * input parameter LAI
	minLeafC = 0.1*treesParams.Al/treesParams.SLA*10000.0;
//base LAI is what you set in the parameter file for lai
	baseLAI = treesParams.Al;
	nsc = getLeafNSC();

		
//when leaves are expanding
//growth respiration is costs are taken from chloroplast sugar and/or starch
//when chloroplast sugar and starch are limiting expansion respiration is
//supported from stored reserves, but leaf growth is slower
//For leaf we convert from the kg N m-2 leaf * lai of m2 m-2 * 10^4 m2 ha-1
	if (delta_lai > 0.0)
	{
//this function adjusts the leaf C/N ratio as a function of available N
//allows C/N to vary from 35-146 (high-low N)
		delta_nsc = RL * delta_lai / treesParams.SLA * 10000.0 * (1.0/0.86);
		CN = 22.0/max(8.0,treesParams.SLA)*35.0 + (1.0-N_avail_rate_plant)*50.0;
		if (delta_nsc > 0.98*nsc)
		{
			delta_nsc = 0.95*nsc;
		}
		GRcost = delta_nsc * 0.14;
		total_N_demand = (delta_nsc-GRcost) / CN;
		if (total_N_demand > 0.99*leafStoredNitrogen[0])
		{
			NfromStorage = 0.99*leafStoredNitrogen[0];
                        delta_lai = NfromStorage * CN * treesParams.SLA / 10000.0;
			delta_nsc = delta_lai / treesParams.SLA * 10000.0 * (1.0/0.86);
			GRcost = delta_nsc * 0.14;
                }
		residual_N_demand = total_N_demand;
		N_neg_demand += N_neg_fract*residual_N_demand;
		N_pos_demand += (1.0-N_neg_fract)*residual_N_demand;

		leafBiomassCarbon[0] += (delta_nsc-GRcost);
		leafNSC[0] -= delta_nsc;
		leafBiomassNitrogen[0] += (delta_nsc-GRcost)/CN;
		leafStoredNitrogen[0] -= (delta_nsc-GRcost)/CN;
	}
//when leaves are senescing and leaf carbon exceeds the leaf bud carbon target
	else if (delta_lai < 0.0 && leafBiomassCarbon[0] > minLeafC)
	{
		delta_nsc = delta_lai / treesParams.SLA * 10000.0;
		CN = leafBiomassCarbon[0]/leafBiomassNitrogen[0];
		leafBiomassCarbon[0] += delta_nsc;
//retain biomass C and N for next budburst
		if (leafBiomassCarbon[0] < minLeafC)
		{
			delta_nsc += (minLeafC-leafBiomassCarbon[0]);
			leafBiomassCarbon[0] = minLeafC;
		}
//0.5 denotes 50% retranslocation 
//apportioned to leaf bud and reserves
		leafBiomassNitrogen[0] += 0.5*delta_nsc/CN;
		leafBiomassNitrogen[0] += 0.25/baseLAI*delta_nsc/CN;
		leafStoredNitrogen[0] -= (0.5-0.25/baseLAI)*delta_nsc/CN;
		leafResidueCarbon[0] -= delta_nsc;
		leafResidueNitrogen[0] -= 0.5*delta_nsc/CN;
	}
}
	
//
//computeLeafNdemand()
//determine how much nitrogen is needed to support leaf growth of delta_lai
//assume that leaf expansion rates require nitrogen to be taken from storage
//any residual nitrogen needed will have to be supplied from belowground
//residual used to compute nitrogen demand
//
void BiogeochemicalCycles::computeLeafNdemand(trees_params treesParams,
						double& delta_lai,
						double N_neg_fract,
						double& N_neg_demand,
						double& N_pos_demand)
{
	double leafCN, total_N_demand, residual_N_demand, NfromStorage;
	double N_avail_rate_plant = plantNstatus[0];

//this function adjusts the leaf C/N ratio as a function of available N
//allows C/N to vary from 35-146 (high-low N)
	leafCN = 22.0/max(8.0,treesParams.SLA)*35.0 + (1.0-N_avail_rate_plant)*50.0;
	total_N_demand = delta_lai / treesParams.SLA * 10000.0 / leafCN;
	NfromStorage = total_N_demand * plantNstatus[0];
	if (NfromStorage > 0.99*leafStoredNitrogen[0])
	{
		NfromStorage = 0.99*leafStoredNitrogen[0];
		delta_lai = NfromStorage * leafCN * treesParams.SLA / 10000.0;
	}
	residual_N_demand = total_N_demand;
	N_neg_demand += N_neg_fract*residual_N_demand;
	N_pos_demand += (1.0-N_neg_fract)*residual_N_demand;
}

//
//updateStemCarbonPools()
//
void BiogeochemicalCycles::updateStemCarbonNitrogenPools(double N_avail_rate_plant,
							 double kratio,
						 	 double nscRatio,
						 	 double rgrowth,
							 double& stemAllocation,
							 double N_neg_fract,
							 double& N_neg_demand,
							 double& N_pos_demand,
						 	 double leafCfraction,
							 double t_canopy,
						 	 double lai,
						 	 double lifeSpan,
						 	 double SLA)
{
	double liveStemIncrement, deadStemIncrement, residueIncrement;
	double tgrowth, CN, stemCincrement;

//increase the rate of live stem death linearly with temperature above 5 C until root warms to 25 C
//reduce root growth at temperatures higher than 25 C
	tgrowth = 0.0;
	if (t_canopy > 5.0)
	{
		tgrowth = (t_canopy-5.0)/20.0;
	}
//Live stem wood mortality
	CN = liveStemCarbon[0]/liveStemNitrogen[0];
	CN = max(CN, 20.0);
	liveStemIncrement = tgrowth*liveStemCarbon[0]/(5.0*lifeSpan);
	liveStemCarbon[0] -= liveStemIncrement;
	liveStemNitrogen[0] -= liveStemIncrement/CN;
	
//Dead stem wood
	deadStemIncrement = liveStemIncrement;
	deadStemCarbon[0] += deadStemIncrement;
	deadStemNitrogen[0] += deadStemIncrement/CN;

//Stem wood residue
	CN = deadStemCarbon[0]/deadStemNitrogen[0];
	CN = max(CN, 20.0);
	residueIncrement = deadStemCarbon[0]*0.01/365.25/48.0;
	stemResidueCarbon[0] += residueIncrement;
	stemResidueNitrogen[0] += residueIncrement/CN;
	deadStemCarbon[0] -= residueIncrement;
	deadStemNitrogen[0] -= residueIncrement/CN;

//increase the rate of growth linearly with temperature above 5 C until root warms to 25 C
//reduce root growth at temperatures higher than 25 C
	tgrowth = 0.0;
	if (t_canopy > 5.0)
	{
		tgrowth = (t_canopy-5.0)/20.0;
		if (tgrowth > 1.0)
		{
			tgrowth = 1.0/tgrowth;
		}
	}
//New stem growth
//when N availability is less than 50% then halt stem growth
	if (N_avail_rate_plant > 0.5 && nscRatio > 1.0)
	{
		stemAllocation = 0.2 * (N_avail_rate_plant-0.5)/0.5 * kratio * tgrowth;
		//stemAllocation = 0.2 * (N_avail_rate_plant-0.5)/0.5 * tgrowth;
	}
	else
	{
		stemAllocation = 0.0;
	}
	stemCincrement = stemAllocation*rgrowth;

//check to make sure there is enough NSC for stem growth
	if (stemCincrement > 0.99*stemNSC[0])
	{
		stemCincrement = 0.99*stemNSC[0];
	}
	stemNSC[0] -= stemCincrement;

	CN = 50.0 + (1.0-N_avail_rate_plant)*50.0;
	N_neg_demand += N_neg_fract*stemCincrement*0.86/CN;
        N_pos_demand += (1.0-N_neg_fract)*stemCincrement*0.86/CN;

	liveStemCarbon[0] += 0.86*stemCincrement;
	liveStemNitrogen[0] += 0.86*stemCincrement/CN;
	leafStoredNitrogen[0] -= 0.86*stemCincrement/CN;
	if (liveStemCarbon[0] < 0.01)
	{
		liveStemCarbon[0] = 0.01;
		liveStemNitrogen[0] = liveStemCarbon[0]/CN;
	}
	if (rgrowth > 0.0)
        {       
                stemAllocation = stemCincrement/rgrowth;
        }
	else
	{
		stemAllocation = 0.0;
	}
}

//
//updateRootCarbonNitrogenPools()
//set root allocation when N is not limiting to decline with root size;
//    start with 0.19, decrease by 0.02 per increment
//set root allocation when N is limiting to shift more C to finer roots
//    start with an average of 0.19 and 0.5, then 0.17 and 0.25, 0.15 and 0.125, 0.13 and 0.125
//    for 10 root classes this sums to 1.0
//DSM - July 2015
//
void BiogeochemicalCycles::updateRootCarbonNitrogenPools(trees_params treesParams,
						 	 double* tempSoil,
						 	 double rgrowth,
							 double N_neg_fract,
							 double& N_neg_demand,
							 double& N_pos_demand,
							 double* kratio_vector,
						 	 double stressedLeafLifeSpan,
						 	 double unstressedLeafLifeSpan,
						 	 double N_avail_rate_plant,
						 	 double leafCfraction,
						 	 double stemAllocation)
{
	double fineRootLow, fineRootHigh, rootLifeSpan, refLifeSpan;
	double referenceSLA, rootScalar, SLAscalar, rootAllocation, rootCincrement, rootCdecrement, residual;
	double CN, tgrowth, kratio_sum, root_relative_growth;

//assumption here is that plants with lower SLA have higher C:N ratios in roots
	referenceSLA = 22.0;
        SLAscalar = referenceSLA / treesParams.SLA;
	if (SLAscalar > 3.0)
	{
		SLAscalar = 3.0;
	}
	else if (SLAscalar < 1.0)
	{
		SLAscalar = 1.0;
	}
	refLifeSpan = treesParams.minRootLifespan*365.25*48.0; //units of 30 min

//sum of dimensionless root K values, each given by K / Ksat
	kratio_sum = 0.0;
        for (int j = 0; j < nRoots; j++)
        {
		if (treesParams.drlat[j+3] < 0.00001)
		{
			kratio_vector[j] = 0.0;
		}
		kratio_sum += kratio_vector[j];
	}
//
//Cycle through all soil-root layers
//
        for (int j = 0; j < nRoots; j++)
        {
//allow for increased C allocation to finer roots at low N
//note that 0.19 + 0.17 +...0.01 = 1.0
//   and 0.5 + 0.25 +...0.00097656 ~= 1.0
        	fineRootLow = 0.19; //allocation to fine roots at optimal available N
        	fineRootHigh = 0.5;
//adjust root lifespan upward as the finest root C:N ratio increases
//root lifespan will also increase with diameter
//differences will be captured between root depths based on C:N
//e.g., McCormack et al 2012. Predicting fine root lifespan from plant functional traits in
//		temperate trees. New Phytologist, 195, 823-831.
		CN = fineRootBiomassCarbon[j][0] / fineRootBiomassNitrogen[j][0];
		CN = max(CN, 20.0);
                //rootLifeSpan = (referenceSLA/treesParams.SLA*refLifeSpan);
                rootLifeSpan = refLifeSpan*CN/20.0;
		assert(rootLifeSpan > 0.000001);
//increase the rate of growth linearly with root temperature above 5 C until root warms to 25 C
//reduce root growth at temperatures higher than 25 C
		tgrowth = 0.0;
		if (tempSoil[j] > 5.0)
		{
			tgrowth = (tempSoil[j]-5.0)/20.0;
			if (tgrowth > 1.0)
			{
				tgrowth = 1.0;
			}
		}
//compute relative allocation of carbon to the respective root-soil layers
		root_relative_growth = kratio_vector[j] / (kratio_sum+0.0000001);
//
//fine root dynamics
//
		rootScalar = 1.0;
                for (int k = 0; k < nFineRootOrders; k++)
                {
//compute root mortality
			CN = fineRootBiomassCarbon[j][k] / fineRootBiomassNitrogen[j][k];
			CN = max(CN, 20.0);
			rootCdecrement = tgrowth*fineRootBiomassCarbon[j][k]/(rootLifeSpan+0.000001);
                        fineRootBiomassCarbon[j][k] -= rootCdecrement;
			fineRootBiomassNitrogen[j][k] -= rootCdecrement/CN;
//change in residue for updating rhizosphere additions (ADD)
			dCrootResidue[j][k] = rootCdecrement;
			dNrootResidue[j][k] = rootCdecrement/CN;
//add dead roots to residue, which will be released over time to the litter C and N
			rootResidueCarbon[j][k] += rootCdecrement;
			rootResidueNitrogen[j][k] += rootCdecrement/CN;

//compute new root growth
//rootAllocation is a weighting of optimal N and stress N allocation fractions
                	rootAllocation = N_avail_rate_plant*fineRootLow;
                        rootAllocation += (1.0-N_avail_rate_plant)*fineRootHigh;
//rootAllocation is further weighted based on hydraulic stress
			rootAllocation *= root_relative_growth;
//convert root NSC and N to root biomass
//CN ratio is increased when the plant has low N reserves
			CN = 20.0*rootScalar*SLAscalar + (1.0-N_avail_rate_plant)*20.0*rootScalar*SLAscalar;
			if (CN > 200.0)
			{
				CN = 200.0;
			}
			rootCincrement = (1.0-leafCfraction-stemAllocation)*rgrowth*tgrowth*rootAllocation;
if (isnan(rootCincrement))
{
cout << leafCfraction << '\t' << stemAllocation << '\t' << rgrowth << '\t' << tgrowth << '\t' << rootAllocation << endl;
exit(1);
}

			if (rootCincrement > 0.99*rootNSC[j][k])
			{
				rootCincrement = 0.99*rootNSC[j][k];
			}
			if (rootCincrement/CN > 0.99 * rootMineralNitrogen[j][k])
			{
				residual = rootCincrement - 0.99 * CN * rootMineralNitrogen[j][k];
				residual *= min(1.0,leafStoredNitrogen[0]/leafBiomassNitrogen[0]*(1.0-leafCfraction));
				rootCincrement = 0.99 * CN * rootMineralNitrogen[j][k] + residual;
				rootCincrement /= 0.86;
				leafStoredNitrogen[0] -= residual/CN;
				rootMineralNitrogen[j][k] *= 0.01;
			}
			else
			{
				rootMineralNitrogen[j][k] -= rootCincrement/CN;
			}

			rootNSC[j][k] -= rootCincrement;
                        fineRootBiomassCarbon[j][k] += 0.86*rootCincrement;
                        fineRootBiomassNitrogen[j][k] += 0.86*rootCincrement/CN;
			N_neg_demand += 0.86*N_neg_fract*rootCincrement/CN;
			N_pos_demand += 0.86*(1.0-N_neg_fract)*rootCincrement/CN;
//increase root lifespan at sqrt(2) with each diameter doubling
                        rootLifeSpan *= treesParams.rootDiamMultiplier; 
                        fineRootLow -= 0.02;
                        fineRootHigh *= 0.5;
			rootScalar *= 1.25;
                }
//
//coarse root dynamics
//
                for (int k = nFineRootOrders; k < nRootOrders; k++)
                {
			CN = coarseRootBiomassCarbon[j][k]/coarseRootBiomassNitrogen[j][k];
			CN = max(CN, 20.0);
//kill off roots
			rootCdecrement = tgrowth*coarseRootBiomassCarbon[j][k]/(rootLifeSpan+0.000001);
                        coarseRootBiomassCarbon[j][k] -= rootCdecrement;
                        coarseRootBiomassNitrogen[j][k] -= rootCdecrement/CN;
//change in residue for updating rhizosphere additions (ADD)
			dCrootResidue[j][k] = rootCdecrement;
			dNrootResidue[j][k] = rootCdecrement/CN;
//add dead roots to residue
			rootResidueCarbon[j][k] += rootCdecrement;
			rootResidueNitrogen[j][k] += rootCdecrement/CN;

//compute new root growth
//rootAllocation is a weighting of optimal N and stress N allocation fractions
                        rootAllocation = N_avail_rate_plant*fineRootLow;
                        rootAllocation += (1.0-N_avail_rate_plant)*fineRootHigh;
//rootAllocation is further weighted based on hydraulic stress
			rootAllocation *= root_relative_growth;
//convert root NSC and N to root biomass
//CN ratio is increased when the plant has low N reserves
			CN = 48.0*SLAscalar + (1.0-N_avail_rate_plant)*48.0*SLAscalar;
			if (CN > 300.0)
			{
				CN = 300.0;
			}
			rootCincrement = (1.0-leafCfraction-stemAllocation)*rgrowth*tgrowth*rootAllocation;
			if (rootCincrement > 0.99*rootNSC[j][k])
			{
				rootCincrement = 0.99*rootNSC[j][k];
			}
			if (rootCincrement/CN > 0.99 * rootMineralNitrogen[j][k])
			{
				rootCincrement = 0.99 * CN * rootMineralNitrogen[j][k];
				rootCincrement /= 0.86;
			}
			rootNSC[j][k] -= rootCincrement;
                        coarseRootBiomassCarbon[j][k] += 0.86*rootCincrement;
                        coarseRootBiomassNitrogen[j][k] += 0.86*rootCincrement/CN;
			rootMineralNitrogen[j][k] -= 0.86*rootCincrement/CN;
			N_neg_demand += 0.86*N_neg_fract*rootCincrement/CN;
			N_pos_demand += 0.86*(1.0-N_neg_fract)*rootCincrement/CN;
//roots double in lifespan with diameter doubling
                        rootLifeSpan *= treesParams.rootDiamMultiplier; 
                        fineRootLow -= 0.02;
                        fineRootHigh *= 0.5;
                }
	}
}

//
//The following is for root maintenance respiration
//   output units of kgC 30min-1
//
double BiogeochemicalCycles::root_respiration_rate(double resp_coef_root,
                                		   double resp_coefficient,
                                		   double t_soil,
                                		   double Croot,
                                		   double transport_rate)
{
        double rate;

        rate = transport_rate * resp_coef_root * exp(resp_coefficient*t_soil) * Croot;

        if (rate < 0.0)
        {
                rate = 0.0;
        }
        return (rate);
}

//
//The following is for stem maintenance respiration
//   output units of kgC 30min-1
//
double BiogeochemicalCycles::stem_respiration_rate(double resp_coef_stem,
                                		   double resp_coefficient,
                                		   double t_canopy,
                                		   double Cstem)
{
        double rate;

        rate = resp_coef_stem * exp(resp_coefficient*t_canopy) * exp(0.67*log(Cstem*10.0))/10.0;

        if (rate < 0.0)
        {
                rate = 0.0;
        }
        return (rate);
}

//
//The following is for leaf maintenance respiration
//   output units of kgC 30min-1
//
double BiogeochemicalCycles::leaf_respiration_rate(double resp_coef_leaf,
                                		   double resp_coefficient,
                                		   double tLeaf,
                                		   double lai,
                                		   double SLA)
{
        double rate;

        rate = resp_coef_leaf * exp(resp_coefficient*tLeaf) * lai / SLA;

        if (rate < 0.0)
        {
                rate = 0.0;
        }
        return (rate);
}

//
//getRhizosphereN()
//
void BiogeochemicalCycles::getRhizosphereN(double& N_neg,
				           double& N_pos)
{
	N_neg = N_pos = 0.0;
        for (int j = 0; j < nRoots; j++)
        {
                for (int k = 0; k < nRootOrders; k++)
                {
                        N_neg += getRhizosphereN_neg(j, k);
                        N_pos += getRhizosphereN_pos(j, k);
                }
        }
}
void BiogeochemicalCycles::getRhizosphereN(double& N_neg,
                                           double& N_pos,
					   int j)
{
        N_neg = N_pos = 0.0;
        for (int k = 0; k < nRootOrders; k++)
        {
                N_neg += getRhizosphereN_neg(j, k);
                N_pos += getRhizosphereN_pos(j, k);
        }
}

//
//getRhizosphereN_neg()
//
double BiogeochemicalCycles::getRhizosphereN_neg(int j,
                                        	 int k)
{
        double Nr;
        Nr = rhizosphereNitrateNitrogen[j][k];
        return(Nr);
}

//
//getRhizosphereN_pos()
//
double BiogeochemicalCycles::getRhizosphereN_pos(int j,
                                        	 int k)
{
        double Nr;
        Nr = rhizosphereAmmoniumNitrogen[j][k];
        return(Nr);
}

//
//getPlantN()
//
double BiogeochemicalCycles::getPlantN()
{
        double plantN = leafStoredNitrogen[0];
        for (int j = 0; j < nRoots; j++)
        {
                for (int k = 0; k < nRootOrders; k++)
                {
                        plantN += rootMineralNitrogen[j][k];
                }
        }
	return(plantN);
}

//
//getRootN()
//
double BiogeochemicalCycles::getRootN()
{
        double rootN = 0.0;
        for (int j = 0; j < nRoots; j++)
        {
                for (int k = 0; k < nRootOrders; k++)
                {
                        rootN += rootMineralNitrogen[j][k];
                }
        }
	return(rootN);
}
double BiogeochemicalCycles::getRootN(int j)
{
        double rootN = 0.0;
        for (int k = 0; k < nRootOrders; k++)
        {
        	rootN += rootMineralNitrogen[j][k];
        }
        return(rootN);
}

//
//getFineRootN()
//
double BiogeochemicalCycles::getFineRootN()
{
        double rootN = 0.0;
        for (int j = 0; j < nRoots; j++)
        {
                for (int k = 0; k < nFineRootOrders; k++)
                {
                        rootN += rootMineralNitrogen[j][k];
                }
        }
	return(rootN);
}
double BiogeochemicalCycles::getFineRootN(int j)
{
        double rootN = 0.0;
        for (int k = 0; k < nFineRootOrders; k++)
        {
        	rootN += rootMineralNitrogen[j][k];
        }
        return(rootN);
}

//
//getRootBiomassN()
//
double BiogeochemicalCycles::getRootBiomassN()
{
        double rootN = 0.0;
        for (int j = 0; j < nRoots; j++)
        {
                for (int k = 0; k < nRootOrders; k++)
                {
                        rootN += fineRootBiomassNitrogen[j][k];
                        rootN += coarseRootBiomassNitrogen[j][k];
                }
        }
	return(rootN);
}
double BiogeochemicalCycles::getRootBiomassN(int j)
{
        double rootN = 0.0;
        for (int k = 0; k < nRootOrders; k++)
        {
        	rootN += fineRootBiomassNitrogen[j][k];
                rootN += coarseRootBiomassNitrogen[j][k];
        }
        return(rootN);
}

//
//getFineRootBiomassN()
//
double BiogeochemicalCycles::getFineRootBiomassN()
{
        double rootN = 0.0;
        for (int j = 0; j < nRoots; j++)
        {
                for (int k = 0; k < nFineRootOrders; k++)
                {
                        rootN += fineRootBiomassNitrogen[j][k];
                        rootN += coarseRootBiomassNitrogen[j][k];
                }
        }
	return(rootN);
}
double BiogeochemicalCycles::getFineRootBiomassN(int j)
{
        double rootN = 0.0;
        for (int k = 0; k < nFineRootOrders; k++)
        {
        	rootN += fineRootBiomassNitrogen[j][k];
                rootN += coarseRootBiomassNitrogen[j][k];
        }
        return(rootN);
}

//
//computeRootNitrogenUptake()
//based on Porporato et al 2003 Advances in Water Resources
//nitrogen uptake from the rhizosphere volume by the root
//passive N uptake is proportional to rhizosphere flux into root
//active N uptake is proportion to ion gradient, for now assuming zero concentration inside root
//rhizopshere flux is already a proportion of canopy transpiration
//proportional to root area / total root area
// Pools are in kg ha-1
// Convert to g m-2 -> Pools * 1000 / 10000
// Water flux is m3 30 min-1
// Root area is per unit ground area
//
void BiogeochemicalCycles::computeRootNitrogenUptake(double UP_neg[][10],
						     double UP_pos[][10],
						     trees_params treesParams,
						     double* thetaSoil,
						     double* Rflux,
						     double ratot,
						     double field_capacity,
						     double DEM_neg,
						     double DEM_pos)
{
	double UP_neg_p, UP_neg_a, UP_pos_p, UP_pos_a, UPp, UPa, DEM_neg_A, DEM_pos_A;
	double F, DEMfract, rootLength, rhizVolume, rhizWaterVolume, kuN, concentration;
	double a, rhizWidth, rootRadius, fluxM, rootFract, absorbFract, totRoot;
	double rootC, totFineRootC, kC, kN, CostUP;
	double DEM_tot, passiveUptake, activeUptake, uptakeRatio;
	double porosity, leafToRootCratio;
	double passiveNitrogenDemand;
	double referenceSLA = 22.0;
        double refLifeSpan = treesParams.minRootLifespan; //years

	porosity = treesParams.porosity;

//this accounts for a gene that controls whether a plant needs to store N from passive uptake
//downregulates uptake as need declines
//start with a linear function
//try a sigmoid function
	//passiveNitrogenDemand = 1.0+ 9.0/(1.0+exp(10.0*plantNstatus[0]-5.0));
	passiveNitrogenDemand = 1.0;

//S.A. Barber, 1995. Soil Nutrient Bioavailability: A mechanistic approach. i
//   John Wiley & Sons, New York
//   Table 4.5 on page 93
//   F = 5 x 10-7 cm2 s-1 * 1800 s 30min-1 * 0.0001 m2 cm-2 = 9 x 10-8 m2 30min-1
//   	F = 9.0e-8; // m2 30min-1;
//
//Jackson 1996, F = 2 x 10-6 cm2 s-1 * 1800 s 30 min-1 * 0.0001 m2 cm-2 = 3.6 x 10-7 m2 30min-1
   	F = 3.6e-7; // m2 30min-1;
//Adjust for emphasis on first-order roots
	F *= 10.0;
//Adjust the root affinity for N uptake as a function of plant N status
	//F /= pow(plantNstatus[0]+0.0000001,2.0);
//proportion of leaf+root biomass in leaf
	totFineRootC = getFineRootCarbon();
	leafToRootCratio = leafBiomassCarbon[0] / (leafBiomassCarbon[0]+totFineRootC);

//rhizosphere width
	rhizWidth = treesParams.rhizosphere_width/1000.0;

//cout << getRootN() << "\t";

//For all soil-root zones
	UP_neg_p = UP_neg_a = UP_pos_p = UP_pos_a = 0.0;
	for (int j = 0; j < nRoots; j++)
	{
		totRoot = 0.0;
		for (int k = 0; k < nRootOrders; k++)
		{
			totRoot += rootArea[j][k];
		}


//convert rhizopshere flux from mmol m-2 s-1 units to m 30min-1 (note: actually m3 30 min-1 m-2)
		fluxM = 0.001*1800.0*EcConversion_StoT(Rflux[j], treesParams.lai);

//either assume nutrients are not allowed to be moved from root to rhizosphere
//or allow rhizodeposition of N
//pick your poison, commenting out this if{} will allow rhizodeposition
		if (fluxM < 0.0)
		{
			fluxM = 0.0;
		}
//For all rhizosphere zones
		//rootRadius = 0.000125;
		//rootRadius = 0.5 *treesParams.maxRootDiam * pow(0.5, 9);
		rootRadius = 0.5 * 0.000125;
//assume 100% of finest root can absorb N, 50% of second finest, etc
		absorbFract = 1.0; 
		for (int k = 0; k < nRootOrders; k++)
		{
			rootFract = absorbFract*rootArea[j][k]/totRoot;
//rhizosphere water volume (m3 m-2)
			rhizVolume = computeRhizosphereVolume(j, k, porosity, rootRadius, rhizWidth);
			rootLength = rootArea[j][k]/(2.0*M_PI*rootRadius);
/*
			if (rhizVolume < 0.000001)
			{
				rhizVolume = 0.000001;
			}
			if (rootLength < 0.001)
			{
				rootLength = 0.001;
			}
*/

//how much liquid water is in the rhizosphere -- this is correct only when there is zero water uptake
			rhizWaterVolume = rhizVolume*thetaSoil[j]/porosity;
			if (rhizWaterVolume < 0.00001)
			{
				rhizWaterVolume = 0.00001;
			}
//
//NITRATE NITROGEN
//
//passive nitrate uptake
//conc = g m-2 / (m3 m-2) = g m-3 water
			a = 1.0; //for mobile nitrate ion
//N solubility declines at lower soil moisture
			if (thetaSoil[j] < field_capacity)
			{
				a *= sqrt(thetaSoil[j]/field_capacity+0.00000001);
			}
			a *= passiveNitrogenDemand;
			concentration = a * rhizosphereNitrateNitrogen[j][k] * 0.1 / rhizWaterVolume;

//UPp = g m-3 x m3 30min-1 m-2 = g m-2 30min-1
			UPp = concentration * fluxM * rootFract;

//determine carbon cost of passive nitrate uptake
//downregulate uptake if NSC is not available
//From Cannell and Thornley, 2000, Annals of Botany
//Assume a respiratory cost of 1 mol NO3- N requires 0.4 mol C in glucose respired
//  In mass units, cost of NO3- uptake is 0.4 * (12/14) = 0.34g glucose per g NO3- N
//  For NH4+ uptake use 0.17 g glucose C per g NH4+ N taken up
//  These numbers assume use of NSC in form of C6H12O6
//Here is the cost of reducing nitrate to ammonium for assimilation
//  cost of nitrate reduction is 8 mol H (mol N)-1
//  glucose C cost of (8 x 6 x 12/24)/14 = 1.72 kg C kg-1 N reduced
			CostUP = 1.72*UPp*10.0*rhizWaterVolume;
			if (CostUP > leafNSC[0] && CostUP > 0.0)
			{
				UPp *= leafNSC[0]/CostUP;
				leafNSC[0] = 0.0;
			}
			else
			{
				leafNSC[0] -= CostUP;
			}
			UP_neg_p += UPp*10.0*rhizWaterVolume;
			UP_neg[j][k] = UPp;
			rootMineralNitrogen[j][k] += UPp*10.0*rhizWaterVolume;
//active nitrate uptake
//convert area-based N demand to volume-based N demand (g m-3 30min-1)
			DEM_neg_A = DEM_neg * 0.1 / rhizWaterVolume;
			DEM_neg_A /= pow(plantNstatus[0]+0.0001, 2.0);
			DEMfract = DEM_neg_A*rootFract;
//diffusion coefficient F is m2 30min-1 / root length x concentration g m-3 = g m-2 30min-1
//kuN = g m-3 x m2 30min-1 x m2 x m-3 = g m-2 30 min-1
//NEED CONCENTRATION GRADIENT (CONCENTRATION - 0)/HALF RHIZOSPHERE RADIUS
/*
			kuN = 0.0;
			if (plantNstatus[0] < 1.0)
			{
				kuN = concentration/(0.5*rhizWidth) * F * rootArea[j][k] * absorbFract;
			}
*/
			kuN = concentration/(0.5*rhizWidth) * F * rootArea[j][k] * absorbFract;
			UPa = 0.0;
			if (DEMfract < UPp)
			{
				UPa = 0.0;
			}
			else if (kuN < (DEMfract-UPp))
			{
				UPa = kuN;
			}
			else if ((kuN > (DEMfract-UPp)) && (DEMfract-UPp) > 0.0)
			{
				UPa = DEM_neg_A - UPp;
			}

//determine carbon cost of active nitrate uptake
//downregulate uptake if NSC is not available
//From Cannell and Thornley, 2000, Annals of Botany
//Assume a respiratory cost of 1 mol NO3- N requires 0.4 mol C in glucose respired
//  In mass units, cost of NO3- uptake is 0.4 * (12/14) = 0.34g glucose per g NO3- N
//  For NH4+ uptake use 0.17 g glucose C per g NH4+ N taken up
//  These numbers assume use of NSC in form of C6H12O6

			CostUP = 0.34*UPa*10.0*rhizWaterVolume;

//active N uptake is limited by available NSC in the root
//cout << "UPa = " << UPa << endl;
			if (CostUP > rootNSC[j][k] && CostUP > 0.0) 
			{
				UPa *= rootNSC[j][k]/CostUP;
				rootNSC[j][k] = 0.0;
			}
			else
			{
				rootNSC[j][k] -= CostUP;
			}
			UP_neg_a += UPa*10.0*rhizWaterVolume;
			UP_neg[j][k] += UPa;

//Add in the cost of reducing nitrate to ammonium for assimilation
//  cost of nitrate reduction is 8 mol H (mol N)-1
//  glucose C cost of (8 x 6 x 12/24)/14 = 1.72 kg C kg-1 N reduced

			CostUP = 1.72*UPa*10.0*rhizWaterVolume;
			leafNSC[0] -= CostUP;

//move active N uptake to root and leaf in equal proportions
//convert from g m-3 to kg ha-1
			rootMineralNitrogen[j][k] += UPa*10.0*rhizWaterVolume; 

//
//AMMONIUM NITROGEN
//
//passive ammonium uptake
			a = 0.1; //assumes that most ammonium is absorbed in soil matrix
//N solubility declines at lower soil moisture
			if (thetaSoil[j] < field_capacity)
			{
				a *= sqrt(thetaSoil[j]/field_capacity+0.00000001);
			}
			a *= passiveNitrogenDemand;
			concentration = a * rhizosphereAmmoniumNitrogen[j][k] * 0.1 / rhizWaterVolume;
			UPp = concentration * fluxM * rootFract;
			UP_pos_p += UPp*10.0*rhizWaterVolume;
			UP_pos[j][k] = UPp;
			rootMineralNitrogen[j][k] += UPp*10.0*rhizWaterVolume;
//active ammonium uptake
			DEM_pos_A = DEM_pos * 0.1 / rhizWaterVolume;
			DEM_pos_A /= pow(plantNstatus[0]+0.0001, 2.0);
			DEMfract = DEM_pos_A*rootFract;
//F is m2 30min-1 x root length x concentration g m-3 = g 30min-1
/*
			kuN = 0.0;
			if (plantNstatus[0] < 1.0)
			{
				kuN = concentration/(0.5*rhizWidth) * F * rootArea[j][k] * absorbFract;
			}
*/
			kuN = concentration/(0.5*rhizWidth) * F * rootArea[j][k] * absorbFract;
			if (DEMfract < UPp)
			{
				UPa = 0.0;
			}
			else if (kuN < (DEMfract-UPp))
			{
				UPa = kuN;
			}
			else if ((kuN > (DEMfract-UPp)) && (DEMfract-UPp) > 0.0)
			{
				UPa = DEM_pos_A - UPp;
			}

//determine carbon cost of active ammonium uptake
//downregulate uptake if NSC is not available
//  For NH4+ uptake use 0.17 g glucose C per g NH4+ N taken up
//  These numbers assume use of NSC in form of C6H12O6

			CostUP = 0.34*UPa*10.0*rhizWaterVolume;

//active N uptake is limited by available NSC in the root
			if (CostUP > rootNSC[j][k] && CostUP > 0.0) 
			{
				UPa *= rootNSC[j][k]/CostUP;
				rootNSC[j][k] = 0.0;
			}
			else
			{
				rootNSC[j][k] -= CostUP;
			}
			UP_pos_a += UPa*10.0*rhizWaterVolume;
			UP_pos[j][k] += UPa;
//move active N uptake to root
//convert from g m-3 to kg ha-1
			rootMineralNitrogen[j][k] += UPa*10.0*rhizWaterVolume; 
			if (UP_neg[j][k] < 0.0)
			{
				UP_neg[j][k] = 0.0;
			}
			if (UP_pos[j][k] < 0)
			{
				UP_pos[j][k] = 0.0;
			}
//radius doubles with each root size class and fraction of root that takes up N halves
			rootRadius *= treesParams.rootDiamMultiplier;
			absorbFract *= 0.5;
		}
	}
//move passive N uptake to leaf
	passiveUptake = UP_neg_p + UP_pos_p;
	//leafStoredNitrogen[0] += leafToRootCratio*passiveUptake;

//cout << getRootN() << endl;

//compute plant nitrogen status
	DEM_tot = DEM_neg + DEM_pos;
	if (treesParams.usePhenology == true) //perennial plant
	{
		plantNstatus[0] = 1.0*leafToRootCratio*treesParams.leafLifeSpan*leafStoredNitrogen[0]/(leafBiomassNitrogen[0]+0.0001);
	}
	else //annual plant
	{
		plantNstatus[0] = 1.0*leafToRootCratio*leafStoredNitrogen[0]/(leafBiomassNitrogen[0]+0.0001);
	}
	//plantNstatus[0] += (1.0-leafToRootCratio)*2.0*refLifeSpan*referenceSLA/treesParams.SLA*getFineRootN()/(getFineRootBiomassN()+0.0001);
	plantNstatus[0] += 1.0*(1.0-leafToRootCratio)*getFineRootN()/(getFineRootBiomassN()+0.0001);
	if (plantNstatus[0] > 1.0)
	{
		plantNstatus[0] = 1.0;
	}
	else if (plantNstatus[0] < 0.001)
	{
		plantNstatus[0] = 0.001;
	}
}

//
//computeLeaching()
//combine drainage proportion by rhizosphere and N concentration
//drain is in units of m 30 min-1
//
void BiogeochemicalCycles::computeLeaching(double LE_neg[][10],
		     			   double LE_pos[][10],
		     			   trees_params treesParams,
		     			   double* thetaSoil,
					   double* depthSoil,
		     			   double* bypassFlow)
{
	double rootRadius, rhizWidth, rhizVolume, totalRhizVolume[10], rhizWaterVolume, a, concentration;
	double porosity, leachate, rhizVolumeProportion, negLE[10], posLE[10];
	rhizWidth = 0.001*treesParams.rhizosphere_width;
	porosity = treesParams.porosity;
//For all soil-root layers
	for (int j = 0; j < nRoots; j++)
	{
		negLE[j] = 0.0;
		posLE[j] = 0.0;
	}
	for (int j = 0; j < nRoots; j++)
	{
		//rootRadius = 0.000125; //meters
		//rootRadius = 0.5 *treesParams.maxRootDiam * pow(0.5, 9);
		rootRadius = 0.5 * 0.000125;
		totalRhizVolume[j] = computeRhizosphereVolume(j, porosity, rootRadius, rhizWidth);
//For all rhizosphere zones
		for (int k = 0; k < nRootOrders; k++)
		{
//total radius is root radius plus rhizosphere radius
//rhizosphere water volume
			rhizVolume = computeRhizosphereVolume(j, k, porosity, rootRadius, rhizWidth);
			if (rhizVolume < 0.000001)
			{
				rhizVolume = 0.000001;
			}
			rhizVolumeProportion = rhizVolume / totalRhizVolume[j];
			rhizWaterVolume = rhizVolume*thetaSoil[j]/porosity;
//nitrate leaching
//units are g m-3
			a = 1.0; //for mobile nitrate ion
			concentration = 0.1 * a * rhizosphereNitrateNitrogen[j][k]/rhizWaterVolume;
			LE_neg[j][k] = rhizVolumeProportion * bypassFlow[j] * concentration;
			negLE[j] += LE_neg[j][k];
//ammonium leaching
			a = 0.1; //assumes that most ammonium is absorbed in soil matrix
			concentration = 0.1 * a * rhizosphereAmmoniumNitrogen[j][k]/rhizWaterVolume;
			LE_pos[j][k] = rhizVolumeProportion * bypassFlow[j] * concentration;
			posLE[j] += LE_pos[j][k];

			rootRadius *= treesParams.rootDiamMultiplier;
		}
//compute net advection of leachate
//negative results indicate net leachate gain; positive indicates net loss
//this assumes that advected leachate is well mixed
		leachate = 0.0;
		if (j > 0)
		{
			for (int k = 0; k < nRootOrders; k++)
			{
				LE_neg[j][k] -= 0.1*negLE[j-1]*totalRhizVolume[j]/depthSoil[j];
				LE_pos[j][k] -= 0.1*posLE[j-1]*totalRhizVolume[j]/depthSoil[j];
				leachate += LE_neg[j][k] + LE_pos[j][k];
			}
		}
//only net loss or no loss from the top root-soil zone
		else
		{
			for (int k = 0; k < nRootOrders; k++)
			{
				leachate += LE_neg[j][k] + LE_pos[j][k];
			}
		}
//we will output the net leaching in kgN ha-1 to be consistent with reporting
		nitrogenLeaching[j] = 10.0*leachate;
	}
}

//
//compute rhizosphere pore volume
//
double BiogeochemicalCycles::computeRhizosphereVolume(int rootNum,
						      double porosity,
						      double rootRadius,
						      double rhizosphereWidth)
{
	double rhizVolume;

	rhizVolume = 0.0;
	for (int k = 0; k < nRootOrders; k++)
	{
		rhizVolume += computeRhizosphereVolume(rootNum, k, porosity, rootRadius, rhizosphereWidth);
		rootRadius *= 2.0;
	}
	return(rhizVolume);
}
double BiogeochemicalCycles::computeRhizosphereVolume(int rootNum,
						      int rootOrder,
						      double porosity,
						      double rootRadius,
						      double rhizosphereWidth)
{
	double totRadius, rootLength, rhizVolume;

	assert(rootNum >= 0);
	assert(rootNum < nRoots);
	assert(rootOrder >= 0);
	assert(rootOrder < nRootOrders);

	totRadius = rhizosphereWidth + rootRadius;
	rootLength = rootArea[rootNum][rootOrder]/(2.0*M_PI*rootRadius);
	rhizVolume = porosity*rootLength*M_PI*(totRadius*totRadius - rootRadius*rootRadius);
//cout << "rhizVolume = " << rhizVolume << endl;
	return(rhizVolume);
}

//
//updateRhizospherePools()
//litter biomass pools
//  dCldt = ADD + BG - DECl
//  BD = kdCb
//  DEC = phi*fds*kl*Cl
//  dNldt = ADD/CNadd + BD/CNb - DECl/CNl
//humus pools
//  dChdt = rh*DECl - DECh
//  DECh = phi * fds * kh * Cb * Ch
//  dNhdt = rh * DECl/CNh - DECh/CNh
//biomass pool
//  dCbdt = (1-rh-rr)*DECl + (1-rr)*DECh - BD
//  dNbdt = (1-rh*CNl/CNh)*DECl/CNl + DECh/CNh - BD/CNb - PHI
//State variables have units of kg ha-1
//Internally these should be converted to g m-3 rhizosphere volume
//Rhizosphere volume given as m3 m-2
//AtoV: kg ha-1 * (1/10000) ha m-2 * 1000 g kg-1 / (m3 m-2) = g m-3
//
void BiogeochemicalCycles::updateRhizospherePools(trees_params treesParams,
						  double* thetaSoil,
						  double* tempSoil,
						  double* depthSoil,
						  double UP_neg[][10],
						  double UP_pos[][10],
						  double LE_neg[][10],
						  double LE_pos[][10])
{
	double PHI, phi, Cb, Ch, CNh, rr, CNb, Cl, CNl, rh;
	double ADD, CNadd, BD, DECl, DECh, dC, dN, dN2, dN3;
	double MIN, IMM_pos, IMM_neg;
	double NIT, fns, kd, kn;
	double theta, t_soil, depth, cumDepth, rootRadius, CN, bulkDensity;
	double dCldt, dNldt, dChdt, dNhdt, dCbdt, dNbdt, dNposdt, dNnegdt;
	double porosity, rhizVolume, rhizWidth, AtoV, AtoVbulk;
	double rootLifeSpan, refLifeSpan;
	double sumRespiration;
	double cumDist;
	double beta = 0.970;

	bulkDensity = treesParams.BD;
	rhizWidth = 0.001*treesParams.rhizosphere_width; //metres
	porosity = treesParams.porosity;
	cumDepth = 0.0;

//For all soil-root zones
	for (int j = 0; j < nRoots; j++)
	{
		refLifeSpan = treesParams.minRootLifespan*48.0*365.25;
//adjust root lifespan upward as the finest root C:N ratio increases
//root lifespan will also increase with diameter
//differences will be captured between root depths based on C:N
//e.g., McCormack et al 2012. Predicting fine root lifespan from plant functional traits in
//              temperate trees. New Phytologist, 195, 823-831.
                CN = fineRootBiomassCarbon[j][0] / fineRootBiomassNitrogen[j][0];
		CN = max(CN, 20.0);
                //rootLifeSpan = (referenceSLA/treesParams.SLA*refLifeSpan);
                rootLifeSpan = refLifeSpan*CN/20.0;
                assert(rootLifeSpan > 0.000001);
		//rootRadius = 0.000125; //meters
		//rootRadius = 0.5 *treesParams.maxRootDiam * pow(0.5, 9);
		rootRadius = 0.5 * 0.000125; //meters
		theta = thetaSoil[j];
		t_soil = tempSoil[j];
		depth = depthSoil[j];
		cumDepth += depth;

//convert kg ha-1 basis to g m-3 bulk soil pore volume for layer i
		cumDist = 1.0-pow(beta, 100.0*(cumDepth-0.5*depth));
		if (cumDist > 0.95)
		{
			cumDist = 0.95;
		}
//cout << cumDist << endl;
		if (treesParams.usePhenology) //perennial plants
		{
			AtoVbulk = 0.1 / depth * bulkDensity / max(treesParams.ar[j+3], 1.0 - cumDist); 
		}
		else //annual plants
		{
			AtoVbulk = 0.1 / depth * bulkDensity;
		}
//humus C and C/N
//assume well-mixed humus throughout soil layer
		Ch = humusCarbon[j];
		CNh = Ch / humusNitrogen[j];
		Ch *= AtoVbulk; //convert to g m-3 soil volume

		sumRespiration = 0.0;

//For all rhizosphere zones
		for (int k = 0; k < nRootOrders; k++)
		{
			rhizVolume = computeRhizosphereVolume(j, k, porosity, rootRadius, rhizWidth);
			if (rhizVolume < 0.000001)
			{
				rhizVolume = 0.000001;
			}
//Convert mass per area units to concentration by volume
//  kg ha-1 * (1/10000) ha m-2 * 1000 g kg-1 / (m3 m-2) = g m-3
// AtoV = (0.1/depth/area) * (depth/area/rhizVolume);
			AtoV = 0.1/rhizVolume;

//Compute residue input C and C/N to the soil and update residue pools
			if (j == 0) //incorporate aboveground residue and root residue
			{
				if (treesParams.usePhenology == true) //perennial plants
				{
//add leaf residue
					dC = 0.1/48.0*leafResidueCarbon[0]/(730.0/(fns+0.0001))*rhizVolume/depth;
					dN = dC * leafResidueNitrogen[0]/leafResidueCarbon[0];
					leafResidueCarbon[0] -= dC;
					leafResidueNitrogen[0] -= dN;
					ADD = dC;
//add stem residue
					dC = 0.1/48.0*stemResidueCarbon[0]/(36500.0/(fns+0.0001))*rhizVolume/depth;
					dN2 = dC * stemResidueNitrogen[0]/stemResidueCarbon[0];
					stemResidueCarbon[0] -= dC;
					stemResidueNitrogen[0] -= dN2;
					ADD += dC;
				}
				else //annual plants
				{
					ADD = 0.0;
					dN = dN2 = 0.0;
				}
//add root residue
				//dC = rootResidueCarbon[j][k]/(rootLifeSpan/(fns+0.0001));
				dC = rootResidueCarbon[j][k]/rootLifeSpan;
				dN3 = dC * rootResidueNitrogen[j][k]/rootResidueCarbon[j][k];
				rootResidueCarbon[j][k] -= dC;
				rootResidueNitrogen[j][k] -= dN3;
				ADD += dC;
				CNadd = ADD/(dN + dN2 + dN3);
//convert from kg ha-1 to g m-3 rhizosphere volume
				ADD *= AtoV; 

			}
			else //incorporate only root residue
			{
				//ADD = rootResidueCarbon[j][k]/(rootLifeSpan/(fns+0.0001));
				ADD = rootResidueCarbon[j][k]/rootLifeSpan;
				dN = ADD * rootResidueNitrogen[j][k]/rootResidueCarbon[j][k];
				CNadd = ADD / dN;
//cout << "ADD = " << ADD;
//cout << "; CNadd = " << CNadd << endl;
				rootResidueCarbon[j][k] -= ADD;
				rootResidueNitrogen[j][k] -= dN;
				ADD *= AtoV; //convert to g m-3 rhizosphere volume
			}
//microbial C and C/N
			Cb = rhizosphereLiveMicrobialCarbon[j][k];
			CNb = Cb / rhizosphereMicrobialNitrogen[j][k];
			Cb *= AtoV; //convert to g m-3 rhizosphere volume

//litter C and C/N
			Cl = rhizosphereCl[j][k];
			CNl = Cl / rhizosphereNl[j][k];
			Cl *= AtoV; //convert to g m-3 rhizosphere volume

//compute decomposition, mineralization-immobilization, and nitrification
//phi is the coefficient of nitrogen sufficiency for bacteria use (Eq. 23 in Porporato 2003 AWR
//
			phi = computeMineralizationImmobilization(MIN, PHI, IMM_pos, IMM_neg, 
							DECl, DECh, rh, rr, kd, kn, fns, j, k, 
							t_soil, theta, AtoV, AtoVbulk, depth, treesParams);
//amount of dead microbial biomass
//a linear function assumes no predation/competition as the population gets crowded
			BD = kd * Cb * max(1.0,log10(Cb)-2.0);
			//BD = kd * Cb * log10(max(10.0,0.5*Cb));
			//BD = kd * Cb;

//
//update the carbon balance in the fast/litter pool (Eq. 4)
//assuming dead biomass going to labile pool declines as CNl declines towards 20
			dCldt = ADD + max(0.0,1.0-2.0*CNh/CNl)*BD - DECl;
			rhizosphereCl[j][k] += dCldt/AtoV; //convert to kg ha-1

//update the nitrogen balance in the fast/litter pool (Eq. 7)
			dNldt = ADD/CNadd + max(0.0,1.0-2.0*CNh/CNl)*BD/CNb - DECl/CNl;
			rhizosphereNl[j][k] += dNldt/AtoV; //convert to kg ha-1

//Update microbial biomass pools
//carbon balance in the biomass (Eq. 11)
//rr represents the fraction of carbon lost as CO2 in respiration
			dCbdt = (1.0 - rh -rr)*DECl + (1.0 - rr)*DECh - BD;
//cout << "j = " << j << "; dCbdt = " << dCbdt << endl;
			rhizosphereLiveMicrobialCarbon[j][k] += dCbdt/AtoV; //convert to kg ha-1

//balance of nitrogen in the biomass (Eq. 12)
//PHI = MIN - (IMM_pos + IMM_neg);
//does not include rr because that is a carbon pathways only
			dNbdt = (1.0 - rh*CNl/CNh)*DECl/CNl + DECh/CNh - BD/CNb - PHI;
			rhizosphereMicrobialNitrogen[j][k] += dNbdt/AtoV; //convert to kg ha-1
//Update ammonium and nitrate in the rhizosphere
			NIT = fns * kn * Cb * rhizosphereAmmoniumNitrogen[j][k]*AtoV;
			dNposdt = MIN - IMM_pos - NIT - LE_pos[j][k] - UP_pos[j][k];
			rhizosphereAmmoniumNitrogen[j][k] += dNposdt/AtoV; //convert to kg ha-1
			if (rhizosphereAmmoniumNitrogen[j][k] < 0.00000001)
			{
				rhizosphereAmmoniumNitrogen[j][k] = 0.00000001;
			}
			dNnegdt = NIT - IMM_neg - LE_neg[j][k] - UP_neg[j][k];
			rhizosphereNitrateNitrogen[j][k] += dNnegdt/AtoV; //convert to kg ha-1
			if (rhizosphereNitrateNitrogen[j][k] < 0.00000001)
			{
				rhizosphereNitrateNitrogen[j][k] = 0.00000001;
			}
			rhizosphereMineralNitrogen[j][k] += (dNposdt + dNnegdt)/AtoV; //convert to kg ha-1
			if (rhizosphereMineralNitrogen[j][k] < 0.00000001)
			{
				rhizosphereMineralNitrogen[j][k] = 0.00000001;
			}

//balance of carbon in the humus (Eq. 8)
			dChdt = rh * DECl + min(1.0,2.0*CNh/CNl)*BD - DECh;

//nitrogen balance equation (Eq. 10)
			//dNhdt = dChdt/CNh;
			dNhdt = (rh * DECl - DECh)/CNh + min(1.0,2.0*CNh/CNl)*BD/CNb;
			humusCarbon[j] += dChdt/AtoV; //convert to kg ha-1
			humusNitrogen[j] += dNhdt/AtoV; //convert to kg ha-1

			rootRadius *= treesParams.rootDiamMultiplier;
			rootLifeSpan *= treesParams.rootDiamMultiplier;

//Note: CO2 release will be added here as rr * (DECl + DECh)
			sumRespiration += rr * (DECl + DECh) / AtoV;
			
		}
		heterotrophicRespiration[j] = sumRespiration;
	}
}

//
//computeMineralizationImmobilization()
//This function implements Eqs. 18 and 23 in Porporato et al 2003. Ad. Water Res, 26, 45-58.
//These equations regulate the dynamics of composition, mineralization, and immobilization.
//PHI = phi * fd(s)Cb{khCh[1/(C/N)h - (1-rr)/(C/N)b] + klCl[1/(C/N)l - rh/(C/N)h - (1-rh-rr)/(C/N)b}
//phi is a non-dimensional factor to take into account poor nitrogen and immobilization is not sufficient
//fd(s) nondimensional factor describing soil moisture effects on decomosition
// CNb = 8, requires CN of 24, 8 C metabolized for each N and 16 respired as CO2
//CNh varies from 10-12, CNb from 8-12, CNl from 20 to 50
//
double BiogeochemicalCycles::computeMineralizationImmobilization(double& MIN,
						   double& PHI,
						   double& IMM_pos,
						   double& IMM_neg,
						   double& DECl,
						   double& DECh,
						   double& rh, 
						   double& rr, 
						   double& kd, 
						   double& kn, 
						   double& fns,
						   int root_num,
						   int root_order,
						   double t_soil,
						   double theta,
						   double AtoV,
						   double AtoVbulk,
						   double depth,
						   trees_params treesParams)
{
	double phi, phi_num, phi_den, Cb, kh, Ch, CNh, CNb, kl, Cl, CNl;
	double ki_pos, ki_neg, N_pos, N_neg, fds, fds_opt;
	double porosity;

//reaction coefficients, daily basis
	kd = treesParams.kd; //D-1
	kn = treesParams.kn; //m3 D-1 gC-1
	kl = treesParams.kl; //m3 D-1 gC-1
	kh = treesParams.kh; //m3 D-1 gC-1
//convert reaction coefficients to 30min-1 basis, time steps used in TREES
//D-1 -> 30min-1 = 1/48 = 0.02083333
	kd *= 0.02083333; //30-min-1
	kn *= 0.02083333; //m3 30-min-1 gC-1
	kl *= 0.02083333; //m3 30-min-1 gC-1
	kh *= 0.02083333; //m3 30-min-1 gN-1
	ki_pos = ki_neg = 1.0 * 0.02083333; //m3 30min-1 gC-1
	
//compute the dependency of decomposition rate on water
	fns = nitrification_water_function(treesParams.theta_opt, theta);

	assert(root_num >= 0);
	assert(root_num < nRoots);
	assert(root_order >= 0);
	assert(root_order < 10);
//live microbial C and N
	Cb = rhizosphereLiveMicrobialCarbon[root_num][root_order];
	CNb = Cb / rhizosphereMicrobialNitrogen[root_num][root_order];
	Cb *= AtoV;
//humus C and N
	Ch = humusCarbon[root_num];
	CNh = Ch / humusNitrogen[root_num];
	Ch *= AtoVbulk; //assumes well-mixed humus throughout soil layer
//litter C and N
	Cl = rhizosphereCl[root_num][root_order];
	CNl = Cl/rhizosphereNl[root_num][root_order];
	Cl *= AtoV;
//mineral N
	N_pos = rhizosphereAmmoniumNitrogen[root_num][root_order]*AtoV;
	N_neg = rhizosphereNitrateNitrogen[root_num][root_order]*AtoV;

//cout << "Cb = " << Cb << "; Ch = " << Ch << "; Cl = " << Cl << endl;

//kl and kh are rates - see soil respiration rate equations
//rr (0 <= rr <= 1-rh is fraction of organic carbon that goes into CO2 production
//rr varies typically from 0.6 to 0.8, increases with labile carbon 
//rh is the isohumic coefficient in range 0.15 to 0.35 depending on plant residues
//assume rh increases with the litter C/N ratio

	rh = min(0.35, CNh/CNl);
	rr = 0.6;
/*
	if (rh < 0.15)
	{
		rh = 0.15;
	}
	else if (rh > 0.35)
	{
		rh = 0.35;
	}
	rr = 0.7 - 0.10 * rhizosphereLabileCarbon[root_num][root_order]*AtoV / 
				(rhizosphereLabileCarbon[root_num][root_order]*AtoV + Cl);
*/
	if (rr > (1.0-rh))
	{
		rr = 1.0-rh;
	}
//This is Eq. 23 in Porporato
	phi_num = ki_pos*N_pos + ki_neg*N_neg;
	phi_den = kh*Ch*(1.0/CNh - (1.0-rr)/CNb) + kl*Cl*(1.0/CNl - rh/CNh - (1.0-rh-rr)/CNb);
	phi = -phi_num/phi_den;

	if (phi_den > 0.0) //net mineralization is occurring
	{
		phi = 1.0;
	}
	if (phi > 1.0) //may be at maximum immobilization
	{
		phi = 1.0;
	}

//In Porporato et al 2003, fds (decomp rate) depends on soil water content
//This is modified here to incorporate the role of temperature as well
//  and this computation is accomplished by using the DAMM (Davidson et al GCB) model
	fds_opt = DAMM_Cflux(treesParams, treesParams.optimal_soil_T, treesParams.theta_opt, 
								treesParams.porosity, depth, Cl);
	if (fds_opt > 0.0)
	{
		fds = DAMM_Cflux(treesParams, t_soil, theta, treesParams.porosity, depth, Cl);
		fds = fds/fds_opt;
	}
	else
	{
		fds = 0.0;
	}
	if (fds > 1.0)
	{
		fds = 1.0;
	}

	DECl = (phi*fds*kl*Cb)*Cl;
	DECh = (phi*fds*kh*Cb)*Ch;

//This is Eq. 18 in Porporato
//PHI = DECh*(1.0/CNh - (1.0-rr)/CNb) + DECl*(1.0/CNl - rh/CNh - (1.0-rh-rr)/CNb);
	PHI = phi*fds*Cb*(kh*Ch*(1.0/CNh - (1.0-rr)/CNb) + kl*Cl*(1.0/CNl - rh/CNh - (1.0-rh-rr)/CNb));
	if (PHI >= 0.0)
	{
		MIN = PHI;
		IMM_pos = IMM_neg = 0.0;
	}
	else
	{
		MIN = 0.0;
		IMM_pos = ki_pos*N_pos/(ki_pos*N_pos+ki_neg*N_neg) * (-1.0*PHI);
		IMM_neg = ki_neg*N_neg/(ki_pos*N_pos+ki_neg*N_neg) * (-1.0*PHI);
	}
	return(phi);
}


//
//nitrification_water_function()
//define the role of degree of staturation on nitrification rate
//Two parts:
//      (1) linear rise to an optimum (parameter theta_opt)
//      (2) less than linear decrease at soil_theta > theta_opt
//DSM November 12, 2008
//
double BiogeochemicalCycles::nitrification_water_function(double theta_opt,
                             			          double theta_soil)
{
        double ts1, ts2, increasing, decreasing, rate;

        if (theta_opt < 0.01)
	{
                 theta_opt = 0.01;
	}
        if (theta_opt > 1.0)
	{
                theta_opt = 1.0;
	}
        ts1 = theta_soil/theta_opt;
        ts2 = theta_opt/theta_soil;
        increasing = pow(ts1, 2.0);
        decreasing = pow(ts2, 2.0);
        rate = 1.0/(0.5*(increasing+decreasing));
        return rate;
}

//
//slow_temp_decomp_rate()
//
double BiogeochemicalCycles::slow_temp_decomp_rate(double t_soil,
                        			   double optimal_soil_T,
                        			   double slow_mineral_rate_parm)
{
        double rate;

        rate = pow(slow_mineral_rate_parm, (t_soil - optimal_soil_T)/10.0);

        return rate;

}

//
//fast_temp_decomp_rate()
//
double BiogeochemicalCycles::fast_temp_decomp_rate(double t_soil,
                        			   double optimal_soil_T,
                        			   double fast_mineral_rate_parm)
{
        double rate;

        rate = pow(fast_mineral_rate_parm, (t_soil - optimal_soil_T)/10.0);

        return rate;

}

//
//DAMM_Cflux()
//DAMM model for Rhet //
//DAMM model native output is mgC m-2 hr-1
//
double BiogeochemicalCycles::DAMM_Cflux(trees_params treesParams,
                  			double t_soil,
                  			double theta,
                  			double porosity,
                  			double depth,
                  			double Clitter)
{
	double Resp, areaCflux;

	Resp = DAMM_reactionVelocity(treesParams, t_soil, theta, porosity, depth, Clitter);

//areaCflux is now in mgC m-2 hr-1, convert to kgC ha-1 30min-1
    	areaCflux = 10000.0 * depth * Resp;
    	areaCflux *= (1.0/2000.0);  //   gC / m2 / 30min
    	areaCflux *= 10.0;        // to kg C / ha / 30min

    	return areaCflux;
}

//
//DAMM_reactionVelocity()
//DAMM model for Rhet //
//We assume this will yield a reaction velocity for each rhizosphere
//
double BiogeochemicalCycles::DAMM_reactionVelocity(trees_params treesParams,
                                        	   double t_soil,
                                        	   double theta,
                                        	   double porosity,
                                        	   double depth,
                                        	   double Clitter)
{
	double R, EaSx, kMsx, xASx;
        double O2airfrac, Dliq, Dgas, kMO2, Clitter_gm2, Clitter_gm3, Clitter_gcm3, exponent;
        double Sx, O2, MMSx, MMO2, VmaxSx, R_Sx;

	R = 8.314472 * pow(10.0,-3.0);

	EaSx = treesParams.EaSx;
	kMsx = treesParams.kMsx;
	xASx = treesParams.xASx;

        O2airfrac = 0.209; //L O2 per L air
        Dliq = 10.97;
        Dgas = 1.67;
        kMO2 = 0.0452;

// need C content in g/cm3
// Clitter is in units of kg ha-1, assume all fast C is in top 15cm
//UNIT CONVERSION HAS A PROBLEM
        Clitter_gm2 = Clitter * 0.1;
        Clitter_gm3 = Clitter_gm2 / depth;
        Clitter_gcm3 = Clitter_gm3 * pow(10.0,-6.0);
        exponent = 4.0/3.0;
        Sx = Clitter_gcm3 * Dliq * pow(theta, 3.0);
        O2 = Dgas * O2airfrac * pow((porosity - theta),exponent);
        MMSx = Sx / (kMsx+Sx);
	MMO2 = O2/(kMO2+O2);
        //VmaxSx = exp(xASx - EaSx/(R * (t_soil + 273.15)));
        VmaxSx = xASx * exp(-EaSx/(R * (t_soil + 273.15)));
        R_Sx = VmaxSx * MMSx * MMO2;
        return R_Sx;
}

