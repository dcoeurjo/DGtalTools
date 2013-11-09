/**
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU Lesser General Public License as
 *  published by the Free Software Foundation, either version 3 of the
 *  License, or  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 **/
/**
 * @file vol2obj.cpp
 * @ingroup Examples
 * @author David Coeurjolly (\c david.coeurjolly@liris.cnrs.fr)
 *
 * @date 2013/10/13
 *
 *
 * This file is part of the DGtal library.
 */

///////////////////////////////////////////////////////////////////////////////
#include <iostream>
#include "DGtal/base/Common.h"
#include "DGtal/helpers/StdDefs.h"
#include "DGtal/io/readers/GenericReader.h"
#include "DGtal/io/boards/Board3D.h"
#include "DGtal/io/DrawWithDisplay3DModifier.h"
#include "DGtal/io/readers/PointListReader.h"

#include "DGtal/images/ImageSelector.h"

#include "DGtal/topology/LightImplicitDigitalSurface.h"
#include "DGtal/topology/DigitalSurface.h"
#include "DGtal/images/imagesSetsUtils/IntervalForegroundPredicate.h"

#include <boost/program_options/options_description.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/variables_map.hpp>

using namespace std;
using namespace DGtal;
using namespace Z3i;

///////////////////////////////////////////////////////////////////////////////
namespace po = boost::program_options;

int main( int argc, char** argv )
{
  // parse command line ----------------------------------------------
  po::options_description general_opt("Allowed options are: ");
  general_opt.add_options()
    ("help,h", "display this message")
    ("input,i", po::value<std::string>(), "vol file (.vol) , pgm3d (.p3d or .pgm3d, pgm (with 3 dims)) file or sdp (sequence of discrete points)" )
    ("output,o", po::value<std::string>(), "Output OBJ filename" )
    ("thresholdMin,m",  po::value<int>()->default_value(0), "threshold min to define binary shape" )
    ("thresholdMax,M",  po::value<int>()->default_value(255), "threshold max to define binary shape" );

  bool parseOK=true;
  po::variables_map vm;
  try
    {
      po::store(po::parse_command_line(argc, argv, general_opt), vm);
    } catch(const std::exception& ex)
    {
      parseOK=false;
      trace.info()<< "Error checking program options: "<< ex.what()<< endl;
    }
  po::notify(vm);
  if( !parseOK || vm.count("help")||argc<=1)
    {
      std::cout << "Usage: " << argv[0] << " [input-file]\n"
                << "Convert a  volume file into OBJ format\n"
                << general_opt << "\n";
      return 0;
    }

  if(! vm.count("input"))
    {
      trace.error() << " The input filename was defined" << endl;
      return 0;
    }
  if(! vm.count("output"))
    {
      trace.error() << " The output filename was defined" << endl;
      return 0;
    }

  string inputFilename = vm["input"].as<std::string>();
  string outputFilename = vm["output"].as<std::string>();
  int thresholdMin = vm["thresholdMin"].as<int>();
  int thresholdMax = vm["thresholdMax"].as<int>();


  typedef ImageSelector<Domain, unsigned char>::Type Image;
  string extension = inputFilename.substr(inputFilename.find_last_of(".") + 1);
  if(extension!="vol" && extension != "p3d" && extension != "pgm3D" &&
     extension != "pgm3d" && extension != "sdp" && extension != "pgm")
    {
      trace.info() << "File extension not recognized: "<< extension << std::endl;
      return 0;
    }

  if(extension=="vol" || extension=="pgm3d" || extension=="pgm3D")
    {
      unsigned int numDisplayed=0;
      Image image = GenericReader<Image>::import (inputFilename );
      trace.info() << "Image loaded: "<<image<< std::endl;
      Domain domain = image.domain();

      IntervalForegroundPredicate<Image> simplePredicate ( image, thresholdMin, thresholdMax );
      KSpace ks;
      bool space_ok = ks.init ( image.domain().lowerBound(),
                                image.domain().upperBound(), true );
      if ( !space_ok )
        {
          trace.error() << "Error in the Khamisky space construction."<<std::endl;
          return 2;
        }

      Board3D<Z3i::Space, Z3i::KSpace> board(ks);

      typedef SurfelAdjacency<KSpace::dimension> MySurfelAdjacency;
      MySurfelAdjacency surfAdj ( true ); // interior in all directions.

      //Set up digital surface.
      typedef LightImplicitDigitalSurface<KSpace, IntervalForegroundPredicate<Image>  > MyDigitalSurfaceContainer;
      typedef DigitalSurface<MyDigitalSurfaceContainer> MyDigitalSurface;
      SCell bel = Surfaces<KSpace>::findABel ( ks, simplePredicate );

      MyDigitalSurfaceContainer* ptrSurfContainer =
        new MyDigitalSurfaceContainer ( ks, simplePredicate, surfAdj, bel );
      MyDigitalSurface digSurf ( ptrSurfContainer ); // acquired

      for(MyDigitalSurface::ConstIterator it = digSurf.begin(), itend = digSurf.end();
          it != itend; ++it)
        board << *it;


      /*      for(Domain::ConstIterator it = domain.begin(), itend=domain.end(); it!=itend; ++it){
        unsigned char  val= image( (*it) );
        if(val<=thresholdMax && val >=thresholdMin){
          board << *it;
          }*/
      board.saveOBJ(outputFilename);
    }
  else
    if(extension=="sdp")
      {
        Board3D<> board;
        vector<Z3i::Point> vectVoxels = PointListReader<Z3i::Point>::getPointsFromFile(inputFilename);
        for(int i=0;i< vectVoxels.size(); i++)
          board << vectVoxels.at(i);

        board.saveOBJ(outputFilename);
      }


  return 0;
}
