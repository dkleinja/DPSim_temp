#include "DPFieldMap.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>

DPFieldMap::DPFieldMap(TString name, TString fieldMapFile, double strength, double zcenter)
{
    using namespace std;

    //load the ascii field map
    ifstream fin(fieldMapFile.Data());
    if(!fin)
    {
        cerr << "Field map input " << fieldMapFile << " does not exist." << endl;
        throw 1;
    }
    cout << " Initializing magnetic field of " << name << " from field map " << fieldMapFile << endl;

    //Load the range and number of bins in each dimension
    string line;

    getline(fin, line);
    stringstream ssx(line);
    ssx >> nx >> xmin >> xmax;
    x0 = (xmax + xmin)/2.;

    getline(fin, line);
    stringstream ssy(line);
    ssy >> ny >> ymin >> ymax;
    y0 = (ymax + ymin)/2.;

    getline(fin, line);
    stringstream ssz(line);
    ssz >> nz >> zmin >> zmax;
    zmin = zmin + zcenter;
    zmax = zmax + zcenter;
    z0 = (zmax + zmin)/2. + zcenter;

    cout << "  its demensions: " << endl;
    cout << "    X: " << nx << " bins, " << xmin << "cm -- " << xmax << "cm." << endl;
    cout << "    Y: " << ny << " bins, " << ymin << "cm -- " << ymax << "cm." << endl;
    cout << "    Z: " << nz << " bins, " << zmin << "cm -- " << zmax << "cm." << endl;

    //Fill the grid
    for(int i = 0; i < 3; ++i)
    {
        grid[i] = new TH3D(Form("%s_%d", name.Data(), i), Form("%s_%d", name.Data(), i),
                           nx, xmin - 0.5*(xmax - xmin)/(nx-1), xmax + 0.5*(xmax - xmin)/(nx-1),
                           ny, ymin - 0.5*(ymax - ymin)/(ny-1), ymax + 0.5*(ymax - ymin)/(ny-1),
                           nz, zmin - 0.5*(zmax - zmin)/(nz-1), zmax + 0.5*(zmax - zmin)/(nz-1));
    }

    while(getline(fin, line))
    {
        double x, y, z, bx, by, bz;

        stringstream ss(line);
        ss >> x >> y >> z >> bx >> by >> bz;

        z = z + zcenter;

        grid[0]->Fill(x, y, z, bx*strength);
        grid[1]->Fill(x, y, z, by*strength);
        grid[2]->Fill(x, y, z, bz*strength);
    }
}

DPFieldMap::~DPFieldMap()
{
    for(int i = 0; i < 3; ++i) delete grid[i];
}

void DPFieldMap::GetFieldValue(const double Point[3], double* Bfield) const
{
    for(int i = 0; i < 3; ++i) Bfield[i] = 0.;
    if(Point[0] < xmin || Point[0] > xmax || Point[1] < ymin || Point[1] > ymax || Point[2] < zmin || Point[2] > zmax) return;

    for(int i = 0; i < 3; ++i) Bfield[i] = grid[i]->Interpolate(Point[0], Point[1], Point[2]);
    //std::cout << grid[0]->GetName() << " " << Point[0]/cm << " " << Point[1]/cm << " " << Point[2]/cm << std::endl;
    //std::cout << "      " << Bfield[0]/tesla << " " << Bfield[1]/tesla << "  " << Bfield[2]/tesla << std::endl;
}
