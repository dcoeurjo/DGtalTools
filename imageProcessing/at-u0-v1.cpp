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
 * @file at-u0-v1.cpp
 * @ingroup Tools
 * @author Jacques-Olivier Lachaud (\c jacques-olivier.lachaud@univ-savoie.fr )
 * Laboratory of Mathematics (CNRS, UMR 5127), University of Savoie, France
 *
 * @date 2016/10/12
 *
 * A tool file named at-u0-v1.
 *
 * This file is part of the DGtal library.
 */

///////////////////////////////////////////////////////////////////////////////
#include <iostream>

#include <sstream>
#include <string>
#include <functional>
#include <boost/format.hpp>

#include <boost/program_options/options_description.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/variables_map.hpp>


#include "DGtal/base/Common.h"
#include "DGtal/helpers/StdDefs.h"
#include "DGtal/io/readers/GenericReader.h"
#include "DGtal/io/writers/GenericWriter.h"

#include "ATu0v1.h"

/**
@page DocATu0v1 at-u0-v1 

@brief Computes a piecewise smooth approximation of an image, by optimizing the Ambrosio-Tortorelli functional.

@writers Marion Foare, Jacques-Olivier Lachaud

@b Usage: at-u0-v1 -i [input.pgm]

@b Usage: at-u0-v1 -i [input.ppm]

Computes the Ambrosio-Tortorelli reconstruction/segmentation of an input image, either grey-level (.pgm) or color image (.ppm).

\f$ AT_e = \int a.(u-g)^2 + v^2 | \grad u|^2 + le.| \grad v|^2 + (l/4e).(1-v)^2 \f$
 
Discretized as (u 0-form, v 1-form, A vertex-edge bdry, B edge-face bdy)

\f$ E(u,v) = a(u-g)^t (u-g) +  u^t A^t diag(v)^2 A^t u + l e v^t (A A^t + B^t B) v + l/(4e) (1-v)^t (1-v) \f$

\b Allowed \b options \b are:

\code
  -h [ --help ]                         display this message
  -i [ --input ] arg                    the input image PGM or PPM filename.
  -o [ --output ] arg (=AT)             the output image basename.
  -l [ --lambda ] arg                   the parameter lambda.
  -1 [ --lambda-1 ] arg (=0.3125)       the initial parameter lambda (l1).
  -2 [ --lambda-2 ] arg (=5.0000000000000002e-05)
                                        the final parameter lambda (l2).
  -q [ --lambda-ratio ] arg (=1.4142135623730951)
                                        the division ratio for lambda from l1 
                                        to l2.
  -a [ --alpha ] arg (=1)               the parameter alpha.
  -e [ --epsilon ] arg                  the initial and final parameter epsilon
                                        of AT functional at the same time.
  --epsilon-1 arg (=2)                  the initial parameter epsilon.
  --epsilon-2 arg (=0.25)               the final parameter epsilon.
  --epsilon-r arg (=2)                  sets the ratio between two consecutive 
                                        epsilon values of AT functional.
  -n [ --nbiter ] arg (=10)             the maximum number of iterations.
  --snr                                 force computation of SNR.
  --image-snr arg                       the input image without deterioration.
  -v [ --verbose ] arg (=0)             the verbose level (0: silent, 1: less 
                                        silent, etc).
\endcode

@image html resATu0v1-cb2-a1_00000-l1_0000000-u.png "AT alpha=1 lambda=1 on carre noise=0.2"

@b example:

\code
./at-u0-v1 -i ../Images/cerclesTriangle64b02.pgm -o AT -a 0.05 -e 1 --lambda-1 0.1 --lambda-2 0.00001
\endcode

@section at_Doc AT


<center>
<table>
<tr>
<td colspan=5>
epsilon scale space
</td>
</tr>
<tr>
  <td> <img height=100px src="resATu0v1-cb2-e2_0-a0_10000-l0_0062092-u0-v1.png"/> </td>
  <td> <img height=100px src="resATu0v1-cb2-e1_0-a0_10000-l0_0062092-u0-v1.png"/> </td>
  <td> <img height=100px src="resATu0v1-cb2-e0_5-a0_10000-l0_0062092-u0-v1.png"/> </td>
  <td> <img height=100px src="resATu0v1-cb2-e0_25-a0_10000-l0_0062092-u0-v1.png"/> </td>
</tr>
<tr>
      <td align = center rowspan="4"> ./build/at-u0-v1 -i Images/carre2Degradesb02.pgm -o cb2 -a 0.1 --lambda 0.006 --epsilon-1 2.0 --epsilon-2 2.0</td>
</tr>
<tr>
      <td align = center rowspan="4"> ./build/at-u0-v1 -i Images/carre2Degradesb02.pgm -o cb2 -a 0.1 --lambda 0.006 --epsilon-1 2.0 --epsilon-2 1.0</td>
</tr>
<tr>
      <td align = center rowspan="4"> ./build/at-u0-v1 -i Images/carre2Degradesb02.pgm -o cb2 -a 0.1 --lambda 0.006 --epsilon-1 2.0 --epsilon-2 0.5</td>
</tr>
<tr>
      <td align = center rowspan="4"> ./build/at-u0-v1 -i Images/carre2Degradesb02.pgm -o cb2 -a 0.1 --lambda 0.006 --epsilon-1 2.0 --epsilon-2 0.25</td>
</tr>
<tr>
<td colspan=5>
alpha scale space
</td>
</tr>
<tr>
    <td> <img height=200px src="resATu0v1-cb2-a1_00000-l1_0000000-u.png"/> </td>
    <td> <img height=200px src="resATu0v1-cb2-a0_50000-l1_0000000-u.png"/> </td>
    <td> <img height=200px src="resATu0v1-cb2-a0_10000-l1_0000000-u.png"/> </td>
    <td> <img height=200px src="resATu0v1-cb2-a0_05000-l1_0000000-u.png"/> </td>
    <td> <img height=200px src="resATu0v1-cb2-a0_01000-l1_0000000-u.png"/> </td>
</tr>
<tr>
        <td align = center rowspan="5"> ./build/at-u0-v1 -i Images/carre2Degradesb02.pgm -o cb2 -a 1.0 --lambda 1.0 --epsilon-1 2.0 --epsilon-2 0.25</td>
</tr>
<tr>
        <td align = center rowspan="5"> ./build/at-u0-v1 -i Images/carre2Degradesb02.pgm -o cb2 -a 0.5 --lambda 1.0 --epsilon-1 2.0 --epsilon-2 0.25</td>
</tr>
<tr>
        <td align = center rowspan="5"> ./build/at-u0-v1 -i Images/carre2Degradesb02.pgm -o cb2 -a 0.1 --lambda 1.0 --epsilon-1 2.0 --epsilon-2 0.25</td>
</tr>
<tr>
        <td align = center rowspan="5"> ./build/at-u0-v1 -i Images/carre2Degradesb02.pgm -o cb2 -a 0.05 --lambda 1.0 --epsilon-1 2.0 --epsilon-2 0.25</td>
</tr>
<tr>
        <td align = center rowspan="5"> ./build/at-u0-v1 -i Images/carre2Degradesb02.pgm -o cb2 -a 0.01 --lambda 1.0 --epsilon-1 2.0 --epsilon-2 0.25</td>
</tr>
<tr>
<td colspan=5>
lambda scale space (lena)
</td>
</tr>
<tr>
<td> <img height=200px src="resATu0v1-lena-370-b02-a0_48000-l0_2000000-u0-v1.png"/> </td>
<td> <img height=200px src="resATu0v1-lena-370-b02-a0_48000-l0_1000000-u0-v1.png"/> </td>
<td> <img height=200px src="resATu0v1-lena-370-b02-a0_48000-l0_0500000-u0-v1.png"/> </td>
<td> <img height=200px src="resATu0v1-lena-370-b02-a0_48000-l0_0250000-u0-v1.png"/> </td>
<td> <img height=200px src="resATu0v1-lena-370-b02-a0_48000-l0_0125000-u0-v1.png"/> </td>
</tr>
<tr>
        <td align = center rowspan="5"> ./build/at-u0-v1 -i Images/lena-370-b02.ppm -o lena -a 0.48 --lambda-1 0.15 --lambda-2 0.0125 -- lambda-ratio 2.0 --epsilon-1 2.0 --epsilon-2 0.25</td>
</tr>
</table>
</center>

*/


using namespace std;
using namespace DGtal;

int main( int argc, char* argv[] )
{
  using namespace Z2i;

  // parse command line ----------------------------------------------
  namespace po = boost::program_options;
  po::options_description general_opt("Allowed options are: ");
  general_opt.add_options()
    ("help,h", "display this message")
    ("input,i", po::value<string>(), "the input image PPM filename." )
    ("output,o", po::value<string>()->default_value( "AT" ), "the output image basename." )
    ("lambda,l", po::value<double>(), "the parameter lambda." )
    ("lambda-1,1", po::value<double>()->default_value( 0.3125 ), "the initial parameter lambda (l1)." ) // 0.3125
    ("lambda-2,2", po::value<double>()->default_value( 0.00005 ), "the final parameter lambda (l2)." )
    ("lambda-ratio,q", po::value<double>()->default_value( sqrt(2) ), "the division ratio for lambda from l1 to l2." )
    ("alpha,a", po::value<double>()->default_value( 1.0 ), "the parameter alpha." )
    ("epsilon,e", po::value<double>(), "the initial and final parameter epsilon of AT functional at the same time." )
    ("epsilon-1", po::value<double>()->default_value( 2.0 ), "the initial parameter epsilon." )
    ("epsilon-2", po::value<double>()->default_value( 0.25 ), "the final parameter epsilon." )
    ("epsilon-r", po::value<double>()->default_value( 2.0 ), "sets the ratio between two consecutive epsilon values of AT functional." )
    ("nbiter,n", po::value<int>()->default_value( 10 ), "the maximum number of iterations." )
    ("snr", "force computation of SNR." )
    ("image-snr", po::value<string>(), "the input image without deterioration." )
    ("pixel-size,p", po::value<int>()->default_value( 1 ), "the pixel size for outputing images (useful when one wants to see the discontinuities v on top of u)." )
    ("verbose,v", po::value<int>()->default_value( 0 ), "the verbose level (0: silent, 1: less silent, etc)." )
    ;

  bool parseOK=true;
  po::variables_map vm;
  try {
    po::store(po::parse_command_line(argc, argv, general_opt), vm);
  } catch ( const exception& ex ) {
    parseOK = false;
    cerr << "Error checking program options: "<< ex.what()<< endl;
  }
  po::notify(vm);
  if ( ! parseOK || vm.count("help")
                 || !vm.count("input")
                 || (vm.count("snr") && !vm.count("image-snr"))
     )
    {
      cerr << "Usage: " << argv[0] << " -i toto.pgm\n"
       << "Computes the Ambrosio-Tortorelli reconstruction/segmentation of an input image."
       << endl << endl
       << " / "
       << endl
       << " | a.(u-g)^2 + v^2 |grad u|^2 + le.|grad v|^2 + (l/4e).(1-v)^2 "
       << endl
       << " / "
       << endl
       << "Discretized as (u 0-form, v 1-form, A vertex-edge bdry, B edge-face bdy)" << endl
       << "E(u,v) = a(u-g)^t (u-g) +  u^t A^t diag(v)^2 A^t u + l e v^t (A A^t + B^t B) v + l/(4e) (1-v)^t (1-v)" << endl
       << endl
       << general_opt << "\n"
       << "Example: ./at-u0-v1 -i ../Images/cerclesTriangle64b02.pgm -o tmp -a 0.05 -e 1 --lambda-1 0.1 --lambda-2 0.00001 -g"
       << endl;
      return 1;
    }
  string f1  = vm[ "input" ].as<string>();
  string f2  = vm[ "output" ].as<string>();
  double l1  = vm[ "lambda-1" ].as<double>();
  double l2  = vm[ "lambda-2" ].as<double>();
  double lr  = vm[ "lambda-ratio" ].as<double>();
  if ( vm.count( "lambda" ) ) l1 = l2 = vm[ "lambda" ].as<double>();
  if ( l2 > l1 ) l2 = l1;
  if ( lr <= 1.0 ) lr = sqrt(2);
  double a   = vm[ "alpha" ].as<double>();
  double e1  = vm[ "epsilon-1" ].as<double>();
  double e2  = vm[ "epsilon-2" ].as<double>();
  if ( vm.count( "epsilon" ) )
    e1 = e2 =  vm[ "epsilon" ].as<double>();
  double er  = vm[ "epsilon-r" ].as<double>();
  int  verb  = vm[ "verbose" ].as<int>();
  int nbiter = vm[ "nbiter" ].as<int>();
  int pix_sz = vm[ "pixel-size" ].as<int>();
  bool color_image = f1.size() > 4 && f1.compare( f1.size() - 4, 4, ".ppm" ) == 0;
  bool grey_image  = f1.size() > 4 && f1.compare( f1.size() - 4, 4, ".pgm" ) == 0;
  if ( ! color_image && ! grey_image ) 
    {
      trace.error() << "Input image file must be either a PGM (grey-level) or a PPM (color) image with these extensions."
                    << endl;
      return 2;
    }

  KSpace K;
  ATu0v1< KSpace > AT( verb );
  Domain domain;

  typedef ImageContainerBySTLVector<Domain, Color> ColorImage;
  typedef ImageContainerBySTLVector<Domain, unsigned char> GreyLevelImage;
  //---------------------------------------------------------------------------
  if ( color_image ) 
    {
      trace.beginBlock("Reading PPM image");
      ColorImage image = PPMReader<ColorImage>::importPPM( f1 );
      domain = image.domain();
      K.init( domain.lowerBound(), domain.upperBound() - Point::diagonal( 1 ), true );
      AT.init( K );
      AT.addInput( image, [] ( Color c ) -> double { return ((double) c.red())   / 255.0; } );
      AT.addInput( image, [] ( Color c ) -> double { return ((double) c.green()) / 255.0; } );
      AT.addInput( image, [] ( Color c ) -> double { return ((double) c.blue())  / 255.0; } );
      trace.endBlock();
    }
  else if ( grey_image ) 
    {
      trace.beginBlock("Reading PGM image");
      GreyLevelImage image = PGMReader<GreyLevelImage>::importPGM( f1 );
      domain = image.domain();
      K.init( domain.lowerBound(), domain.upperBound() - Point::diagonal( 1 ), true );
      AT.init( K );
      AT.addInput( image, [] (unsigned char c ) { return ((double) c) / 255.0; } );
      trace.endBlock();
    }
  //---------------------------------------------------------------------------
  // Prepare output domain
  Domain out_domain( pix_sz * domain.lowerBound(), 
                     pix_sz * domain.upperBound() + Point::diagonal( pix_sz - 1) );
  //---------------------------------------------------------------------------
  AT.setUFromInput();
  AT.setAlpha( a );
  trace.info() << AT << std::endl;
  double n_v = 0.0;
  double eps = 0.0;
  while ( l1 >= l2 )
    {
      trace.info() << "************ lambda = " << l1 << " **************" << endl;
      AT.setLambda( l1 );
      for ( eps = e1; eps >= e2; eps /= er )
        {
          trace.info() << "  ======= epsilon = " << eps << " ========" << endl;
          AT.setEpsilon( eps );
          int n = 0;
          do {
            trace.progressBar( n, nbiter );
            AT.solveU();
            AT.solveV();
            AT.checkV();
            n_v = AT.computeVariation();
          } while ( ( n_v > 0.0001 ) && ( ++n < nbiter ) );
          trace.progressBar( n, nbiter );
          trace.info() << "[#### last variation = " << n_v << " " << endl;
        }
      if ( grey_image ) 
        {
          if ( verb > 0 ) trace.beginBlock("Writing u[0] as PGM image");
          ostringstream ossU, ossV;
          ossU << boost::format("%s-a%.5f-l%.7f-u.pgm") % f2 % a % l1;
          ossV << boost::format("%s-a%.5f-l%.7f-u-v.pgm") % f2 % a % l1;
          GreyLevelImage image( out_domain );
          functions::dec::dualForm2ToGreyLevelImage
            ( AT.calculus, AT.primal_h0 * AT.getU( 0 ), image, 0.0, 1.0, pix_sz ); 
          PGMWriter<GreyLevelImage>::exportPGM( ossU.str(), image );
          functions::dec::dualForm1ToGreyLevelImage
            ( AT.calculus, AT.primal_h1 * AT.getV(), image, 0.0, 1.0, pix_sz ); 
          PGMWriter<GreyLevelImage>::exportPGM( ossV.str(), image );
          if ( verb > 0 ) trace.endBlock();
        }
      else if ( color_image )
        {
          if ( verb > 0 ) trace.beginBlock("Writing u[0,1,2] as PGM image");
          ostringstream ossU;
          ossU << boost::format("%s-a%.5f-l%.7f-u.ppm") % f2 % a % l1;
          string str_image_u = ossU.str();
          ColorImage image( out_domain );
          functions::dec::threeDualForms2ToRGBColorImage
            ( AT.calculus, 
              AT.primal_h0 * AT.getU( 0 ), AT.primal_h0 * AT.getU( 1 ), AT.primal_h0 * AT.getU( 2 ),
              image, 0.0, 1.0, pix_sz ); 
          PPMWriter<ColorImage, functors::Identity >::exportPPM( str_image_u, image );
          if ( verb > 0 ) trace.endBlock();
        }
      l1 /= lr;
    }
  return 0;
}
