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
 * @file qglViewer.cpp
 * @ingroup Examples
 * @author Bertrand Kerautret (\c kerautre@loria.fr )
 * LORIA (CNRS, UMR 7503), University of Nancy, France
 *
 * @date 2011/01/04
 *
 * An example file named qglViewer.
 *
 * This file is part of the DGtal library.
 */

///////////////////////////////////////////////////////////////////////////////
#include <iostream>
#include <queue>
#include <QImageReader>
#include <QtGui/qapplication.h>
#include "DGtal/io/viewers/Viewer3D.h"
#include "DGtal/io/DrawWithDisplay3DModifier.h"
#include "DGtal/io/Color.h"
#include "DGtal/shapes/Shapes.h"
#include "DGtal/helpers/StdDefs.h"

#include "DGtal/io/readers/VolReader.h"
#include "DGtal/images/imagesSetsUtils/SetFromImage.h"
#include "DGtal/images/imagesSetsUtils/SimpleThresholdForegroundPredicate.h"
#include "DGtal/images/ImageSelector.h"

#include <boost/program_options/options_description.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/variables_map.hpp>


#include "DGtal/geometry/volumes/distance/SeparableMetricHelper.h"
#include "DGtal/geometry/volumes/distance/DistanceTransformation.h"

#include <algorithm>                 // for std::for_each
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/adjacency_list.hpp>

///////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace DGtal;
using namespace Z3i;

namespace po = boost::program_options;

///////////////////////////////////////////////////////////////////////////////

/**
 * Missing parameter error message.
 *
 * @param param
 */
void missingParam ( std::string param )
{
    trace.error() <<" Parameter: "<<param<<" is required..";
    trace.info() <<std::endl;
    exit ( 1 );
}



int main( int argc, char** argv )
{

  // parse command line ----------------------------------------------
  po::options_description general_opt ( "Allowed options are: " );
  general_opt.add_options()
    ( "help,h", "display this message." )
    ( "input,i", po::value<std::string>(), "Input vol file." )
    ( "output,o", po::value<std::string>(), "Ouput filename without any extension" )
    ( "vol,v",  "Export the skeleton as a vol file (notyetimplemented)" )
    ( "graph,g", "Export the skeleton as a weighted graph" )
    ( "viewer", "Display the skeleton using DGtal::Viewer3D" );
 

  po::variables_map vm;
  po::store ( po::parse_command_line ( argc, argv, general_opt ), vm );
  po::notify ( vm );
  if ( vm.count ( "help" ) ||argc<=1 )
    {
      trace.info() << "Illustration of guided homotopic thinning of a vol file with 3D viewer."<<std::endl
                   << std::endl << "Basic usage: "<<std::endl
                   << "\thomotopicThinning3d [options] --input <volFileName>"<<std::endl
                   << general_opt << "\n";
      return 0;
    }

  //Parse options
  if ( ! ( vm.count ( "input" ) ) ) missingParam ( "--input" );
  std::string filename = vm["input"].as<std::string>();
 
  std::string outputfilename  ;
  if  ( vm.count ( "output" ) )
    outputfilename = vm["output"].as<std::string>();
  
  if  ( ( vm.count ( "vol" )  ||  (vm.count("graph"))) && ( ! vm.count("output")))
    missingParam ( "--output (mandatory for either vol or graph export)" );
 
  if  (  vm.count ( "output" )  && (!  (vm.count("graph") ||   vm.count("vol"))))
    missingParam ( "--graph or --vol (mandatory with --output option)" );

  //Starting the process   
  typedef ImageSelector < Z3i::Domain, unsigned char>::Type Image;
  Image image = VolReader<Image>::importVol ( filename );

  trace.beginBlock("DT Computation");
  typedef  DistanceTransformation<Image, 0> DTL2;
  DTL2 dtL2;
  
  DTL2::OutputImage resultL2 = dtL2.compute ( image );
  trace.endBlock();
  trace.info() <<image<<std::endl;

  // Domain cretation from two bounding points.
  Point c( 0, 0, 0 );
  Point p1( -50, -50, -50 );
  Point p2( 50, 50, 50 );
  Domain domain( p1, p2 );
  
  trace.beginBlock("Constructing Set");
  DigitalSet shape_set( domain );
  SetPredicate<DigitalSet> set3dPredicate( shape_set );
  SetFromImage<DigitalSet>::append<Image>(shape_set, image,
                                          0, 255);
  trace.info() << shape_set<<std::endl;
  trace.endBlock();

  trace.beginBlock("Computing skeleton");
  Object26_6 shape( dt26_6, shape_set );
  int nb_simple=0; 
  int layer = 1;
  std::queue<DigitalSet::Iterator> Q;
  do 
    {
      trace.info() << "Layer: "<< layer << std::endl;
      int nb=0;
      DigitalSet & S = shape.pointSet();
 
      for ( DigitalSet::Iterator it = S.begin(); it != S.end(); ++it )
        {
	  trace.progressBar((double)nb, (double)S.size()); 
	  trace.info() << nb<<" "; nb++;
	  if (resultL2( *it ) <= layer*layer)
	    {
	      if ( shape.isSimple( *it ) )
		Q.push( it );
	    }
	}
      nb_simple = 0;
      while ( ! Q.empty() )
        {
          DigitalSet::Iterator it = Q.front();
          Q.pop();
          if ( shape.isSimple( *it ) )
            {
              S.erase( *it );
              ++nb_simple;
            }
        }
      trace.info() << "Nb simple points : "<<nb_simple<<std::endl;
      ++layer;
     }
  while ( nb_simple != 0 );
  trace.endBlock();
  
  trace.info() << "Skeleton --> "<<shape<<std::endl;

  //Graph construction & export
  if (vm.count("graph"))
    {
      trace.beginBlock("Exporting the graph");
      typedef boost::property<boost::vertex_distance_t, double, 
                              boost::property<boost::vertex_name_t, std::string> > VertexProperty;
      
      typedef boost::adjacency_list<boost::setS, boost::vecS, boost::undirectedS
                                     > Graph;
      typedef boost::graph_traits<Graph>::vertex_descriptor Vertex;
      typedef std::pair<int, int> Edge;
      typedef std::set < Edge > Edges;
      std::map < Point, std::size_t > myVertexIndex;
      std::vector < Point > myCoords;
   
      Edges edges;
      
      const int num_nodes = shape.pointSet().size();
      Graph G(num_nodes );
           
      size_t pos = 0;
      for(Object26_6::DigitalSet::ConstIterator it = shape.pointSet().begin(),
            itend = shape.pointSet().end(); it != itend;
          ++it, ++pos)
        {
          myCoords.push_back(*it);
          myVertexIndex[*it] = pos;
        }
      
      pos = 0;
      for(Object26_6::DigitalSet::ConstIterator it = shape.pointSet().begin(),
            itend = shape.pointSet().end(); it != itend; ++it, ++pos)
        {
          Object26_6::SmallObject neig = shape.properNeighborhood( * it );
          for(Object26_6::SmallObject::DigitalSet::ConstIterator itn = neig.pointSet().begin(),
                itnend=neig.pointSet().end(); itn != itnend; ++itn)
            add_edge( pos, myVertexIndex.find(*itn)->second , G);        
        }
      
      //Exporting the skeleton
      std::ofstream out((outputfilename+".graph").c_str());
      
      // get the property map for vertex indices
      typedef boost::property_map<Graph, boost::vertex_index_t>::type IndexMap;
      IndexMap index = get(boost::vertex_index, G);
      
      out << num_nodes<<" "<< boost::num_edges(G)<< std::endl;
      typedef boost::graph_traits<Graph>::vertex_iterator vertex_iter;
      std::pair<vertex_iter, vertex_iter> vp;
      for (vp = boost::vertices(G); vp.first != vp.second; ++vp.first)
        out << index[*vp.first] <<  " "
            << myCoords[index[*vp.first]][0] <<" "
            << myCoords[index[*vp.first]][1] <<" "
            << myCoords[index[*vp.first]][2] <<" "
            << (double) sqrt(resultL2( myCoords[index[*vp.first]] ) ) << std::endl;
      out << std::endl;
      // ...
      
      boost::graph_traits<Graph>::edge_iterator ei, ei_end;
      for (tie(ei, ei_end) = boost::edges(G); ei != ei_end; ++ei)
        out<< index[source(*ei, G)] << " "<< index[target(*ei, G)] << std::endl;
      // ...

      out.close();
      trace.endBlock();
    }
          
        
  // Display by using two different list to manage OpenGL transparency.
  if (vm.count("viewer"))
    {
      QApplication application(argc,argv);
      Viewer3D viewer;
      viewer.setWindowTitle("simpleExample3DViewer");
      viewer.show();  
      
      viewer << SetMode3D( shape_set.className(), "Paving" );
      viewer << CustomColors3D(Color(25,25,255, 255), Color(25,25,255, 255));
      viewer << shape.pointSet() ; 

      viewer << SetMode3D( shape_set.className(), "PavingTransp" );
      viewer << CustomColors3D(Color(250, 0,0, 25), Color(250, 0,0, 5));
      viewer << shape_set;
      
      viewer<< Viewer3D::updateDisplay;
      
      return application.exec();
    }
}
//                                                                           //
///////////////////////////////////////////////////////////////////////////////


