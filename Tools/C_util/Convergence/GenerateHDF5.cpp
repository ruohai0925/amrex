
#include <fstream>
#include <iostream>
#include <string>
// using std::string;
#include <AMReX_ParmParse.H>
#include <AMReX_PlotFileUtil.H>
#include <AMReX_PlotFileUtilHDF5.H>

using namespace amrex;

int
main (int   argc,
      char* argv[])
{
    amrex::Initialize(argc,argv);
    {

        if (argc == 1) {
        //    PrintUsage(argv[0]);
        }

        // plotfile names for the coarse and output
        std::string iFile1;
        std::string outFile="out1";

        // read in parameters from inputs file
        ParmParse pp;

        // coarse MultiFab
        pp.query("infile1", iFile1);
        if (iFile1.empty())
            amrex::Abort("You must specify `infile1'");

        // Read the info from the Header file
        std::string HeaderFile = "";
        if (amrex::FileExists(iFile1+"/Header")) {
            HeaderFile = iFile1 + "/Header";
        }

        // Open the file
        std::ifstream infile(HeaderFile);

        std::string Version = "";
        getline(infile, Version);
        std::cout << "Version " << Version << std::endl;

        int n_data_items = 0;
        infile >> n_data_items;
        // Skip the rest of the line (including the newline character)
        infile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
        std::cout << "n_data_items " << n_data_items << std::endl;

        Vector<std::string> varNames(n_data_items);
        for (int i = 0; i < n_data_items; i++)
        {
            getline(infile, varNames[i]);
            std::cout << "i " << i << " " << varNames[i] << std::endl;            
        }

        int dimension = 2;
        infile >> dimension;
        AMREX_ASSERT(dimension==AMREX_SPACEDIM);
        std::cout << "dimension " << dimension << std::endl;

        Real cumtime = 0.0;
        infile >> cumtime;
        std::cout << "cumtime " << cumtime << std::endl;

        int f_lev = 2;
        infile >> f_lev;
        std::cout << "f_lev " << f_lev << std::endl;

        Vector<Real> problo(dimension), probhi(dimension);
        for (int i = 0; i < dimension; i++)
        {
            infile >> problo[i];
            std::cout << "i " << i << " " << problo[i] << std::endl;             
        }

        for (int i = 0; i < dimension; i++)
        {
            infile >> probhi[i];
            std::cout << "i " << i << " " << probhi[i] << std::endl;             
        }

        // Close the file
        infile.close();

        // single-level for now
        // AMR comes later, where we iterate over each level in isolation

        // check to see whether the user pointed to the plotfile base directory
        // or the data itself
        if (amrex::FileExists(iFile1+"/Level_0/Cell_H")) {
            iFile1 += "/Level_0/Cell";
        }

        // storage for the input coarse and fine MultiFabs
        MultiFab mf_c;

        // read in plotfiles, 'coarse' to MultiFab
        VisMF::Read(mf_c, iFile1);

        if (mf_c.contains_nan()) {
            Abort("First plotfile contains NaN(s)");
        }

        // check number of components
        int ncomp = mf_c.nComp();
        Print() << "ncomp = " << ncomp << std::endl;

        // check nodality
        IntVect c_nodality = mf_c.ixType().toIntVect();
        Print() << "nodality " << c_nodality << std::endl;

        // get coarse boxArray
        BoxArray ba_c = mf_c.boxArray();

        // minimalBox() computes a single box to enclose all the boxes
        // enclosedCells() converts it to a cell-centered Box
        Box bx_c = ba_c.minimalBox().enclosedCells();
        Print() << "Box bx_c = " << bx_c << std::endl;

        // number of cells in the coarse domain
        Print() << "npts in coarse domain = " << bx_c.numPts() << std::endl;

        // grab the distribution map from the coarse MultiFab
        DistributionMapping dm = mf_c.DistributionMap();

        // Here we only wannt output the u, v, and p in the HDF5 file for training.
        int n_out_comp = dimension + 1;
        MultiFab mf_out(ba_c,dm,n_out_comp,0);

        // MultiFab::Copy(mfdst, mfsrc, sc, dc, nc, ng); // Copy from mfsrc to mfdst
        Vector<std::string> varNames_out(n_out_comp);
        for (int i = 0; i < n_data_items; i++)
        {
            if (varNames[i] == "x_velocity") {
                std::cout << "i " << i << " " << varNames[i] << std::endl; 
                MultiFab::Copy(mf_out, mf_c, i, 0, 1, 0);
                varNames_out[0] = varNames[i];
            }
            if (varNames[i] == "y_velocity") {
                std::cout << "i " << i << " " << varNames[i] << std::endl; 
                MultiFab::Copy(mf_out, mf_c, i, 1, 1, 0);
                varNames_out[1] = varNames[i];
            }
            if (varNames[i]=="avg_pressure") {
                std::cout << "i " << i << " " << varNames[i] << std::endl; 
                MultiFab::Copy(mf_out, mf_c, i, n_out_comp-1, 1, 0);
                varNames_out[n_out_comp-1] = varNames[i];
            }
        }

        // // write out the outFile if it was specified at the command line
        if (outFile != "") {

            // define the problem domain
            RealBox real_box;
            real_box = RealBox(&(problo[0]),&(probhi[0]));

            Vector<int> is_periodic(AMREX_SPACEDIM,1);

            // build a geometry object so we can use WriteSingleLevelPlotfile
            Geometry geom(bx_c,&real_box,CoordSys::cartesian,is_periodic.data());

            // // give generic variables names for now
            // Vector<string> varNames(ncomp);
            // for (int i=0; i<ncomp; ++i) {
            //     varNames[i] = std::to_string(i);
            // }

            WriteSingleLevelPlotfileHDF5MultiDset(outFile,mf_out,varNames_out,geom,0.,0);
            // WriteSingleLevelPlotfileHDF5MultiDset("out2",mf_c,varNames,geom,0.,0);
        }

    }
    amrex::Finalize();
}
