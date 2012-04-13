#include "DGtal/io/colormaps/GrayscaleColorMap.h"
#include "DGtal/io/readers/PNMReader.h"


#include "DGtal/images/ImageContainerBySTLVector.h"
#include "DGtal/images/ImageSelector.h"

#include "DGtal/geometry/curves/representation/FreemanChain.h"
#include "DGtal/geometry/helpers/ContourHelper.h"

#include "DGtal/topology/helpers/Surfaces.h"

#include <boost/program_options/options_description.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/variables_map.hpp>

#include <vector>
#include <string>

using namespace DGtal;




///////////////////////////////////////////////////////////////////////////////
namespace po = boost::program_options;



void saveAllContoursAsFc(std::vector< std::vector< Z2i::Point >  >  vectContoursBdryPointels, uint minSize){
  for(unsigned int k=0; k<vectContoursBdryPointels.size(); k++){
    if(vectContoursBdryPointels.at(k).size()>minSize){
      	  FreemanChain<Z2i::Integer> fc (vectContoursBdryPointels.at(k));    
	  cout << fc.x0 << " " << fc.y0   << " " << fc.chain << endl; 
	  
    }
  }
}


void saveSelContoursAsFC(std::vector< std::vector< Z2i::Point >  >  vectContoursBdryPointels, 
			 uint minSize, Z2i::Point refPoint, double selectDistanceMax){

  for(unsigned int k=0; k<vectContoursBdryPointels.size(); k++){
    if(vectContoursBdryPointels.at(k).size()>minSize){
      Z2i::Point ptMean = ContourHelper::getMeanPoint(vectContoursBdryPointels.at(k));
      unsigned int distance = (unsigned int)ceil(sqrt((double)(ptMean[0]-refPoint[0])*(ptMean[0]-refPoint[0])+
						      (ptMean[1]-refPoint[1])*(ptMean[1]-refPoint[1])));
      if(distance<=selectDistanceMax){
	FreemanChain<Z2i::Integer> fc (vectContoursBdryPointels.at(k));    
	cout << fc.x0 << " " << fc.y0   << " " << fc.chain << endl; 
      }      
    }    
  }
}




int main( int argc, char** argv )
{
  // parse command line ----------------------------------------------
  po::options_description general_opt("Allowed options are: ");
  general_opt.add_options()
    ("help,h", "display this message")
    ("image,i", po::value<std::string>(), "image file name")
    ("min,m", po::value<int>(), "min image threshold value (default 128)")
    ("max,M", po::value<int>(), "max image threshold value (default 255)")
    
    ("minSize,s", po::value<int>(), "minSize of the extracted freeman chain (default 0)")
    ("contourSelect,s", po::value<vector <int> >()->multitoken(), 
     "Select contour according reference point and maximal distance:  ex. --contourSelect X Y distanceMax")
    ("thresholdRangeMin,r", po::value<vector <int> >()->multitoken(), 
     "use a range interval as threshold (from min) : --thresholdRangeMin min increment max : for each possible i, it define a digital sets [min, min+((i+1)*increment)] such that min+((i+1)*increment)< max  and extract their boundary. ")
    ("thresholdRangeMax,R", po::value<vector <int> >()->multitoken(), 
     "use a range interval as threshold (from max) : --thresholdRangeMax min increment max : for each possible i, it define a digital sets [ max-((i)*increment), max] such that max-((i)*increment)>min  and extract their boundary. ");
  
  
  
  po::variables_map vm;
  po::store(po::parse_command_line(argc, argv, general_opt), vm);  
  po::notify(vm);    
  if(vm.count("help")||argc<=1)
    {
      trace.info()<< "Extract FreemanChains from thresholded image" <<std::endl << "Basic usage: "<<std::endl
      << "\t image2freeman [options] --image <imageName> -min 128 -max 255 > contours.fc"<<std::endl
      << general_opt << "\n";
      return 0;
    }
  
  
  double minThreshold = 128;
  double maxThreshold = 255;
  unsigned int minSize =0;
  bool select=false;
  bool thresholdRange=vm.count("thresholdRangeMin")||vm.count("thresholdRangeMax");
  Z2i::Point selectCenter;
  unsigned int selectDistanceMax = 0; 
  

  //Parse options
  if (!(vm.count("image"))){
    trace.info() << "Image file name needed"<< endl;
    return 0;
  } 
  
  if(vm.count("min")){
    minThreshold= vm["min"].as<int>();
  } 
  if(vm.count("max")){
    maxThreshold= vm["max"].as<int>();
  } 
  if(vm.count("minSize")){
    minSize = vm["minSize"].as<int>();
  } 
  if(vm.count("contourSelect")){
    select=true;
    vector<int> cntConstraints= vm["contourSelect"].as<vector <int> >();
    if(cntConstraints.size()!=3){
      trace.info() << "Incomplete option \"--contourSelect\""<< endl;
      return 0;
    }
    selectCenter[0]= cntConstraints.at(0);
    selectCenter[1]= cntConstraints.at(1);
    selectDistanceMax= (unsigned int) cntConstraints.at(2);
  }

  int min, max, increment;
  if(! thresholdRange){
    min=(int)minThreshold;
    max= (int)maxThreshold;
    increment =  (int)(maxThreshold- minThreshold);
  }else{
    vector<int> vectRange;
    if ( vm.count("thresholdRangeMax")){
      vectRange= vm["thresholdRangeMax"].as<vector <int> >();
    }else{
      vectRange= vm["thresholdRangeMin"].as<vector <int> >();
    }
    if(vectRange.size()!=3){
      trace.info() << "Incomplete option \"--thresholdRange\""<< endl;
      return 0;
    }
    min=vectRange.at(0);
    increment=vectRange.at(1);
    max = vectRange.at(2);
    minThreshold=min;
    maxThreshold=max;
  }

 

  
  typedef ImageSelector < Z2i::Domain, unsigned char>::Type Image;
  typedef IntervalThresholder<Image::Value> Binarizer; 
  string imageFileName = vm["image"].as<std::string>();
  Image image = PNMReader<Image>::importPGM( imageFileName ); 
  
  Z2i::KSpace ks;
  if(! ks.init( image.domain().lowerBound(), 
		image.domain().upperBound(), true )){
    trace.error() << "Problem in KSpace initialisation"<< endl;
  }
  
  
  if (!thresholdRange){
    Binarizer b(min, max); 
    PointFunctorPredicate<Image,Binarizer> predicate(image, b); 
    trace.info() << "DGtal contour extraction from thresholds ["<<  min << "," << max << "]" ;
    SurfelAdjacency<2> sAdj( true );
    std::vector< std::vector< Z2i::Point >  >  vectContoursBdryPointels;
    Surfaces<Z2i::KSpace>::extractAllPointContours4C( vectContoursBdryPointels,
						      ks, predicate, sAdj );  
    if(select){
      saveSelContoursAsFC(vectContoursBdryPointels,  minSize, selectCenter,  selectDistanceMax);
    }else{
      saveAllContoursAsFc(vectContoursBdryPointels,  minSize); 
    }
  }else{
    for(int i=0; minThreshold+i*increment< maxThreshold; i++){
      if(vm.count("thresholdRangeMin")){
	min = (int)(minThreshold+(i)*increment);
      }
      if(vm.count("thresholdRangeMax")){
	max = (int)(maxThreshold-(i)*increment);
      }
      Binarizer b(min, max); 
      PointFunctorPredicate<Image,Binarizer> predicate(image, b); 
      
      trace.info() << "DGtal contour extraction from thresholds ["<<  min << "," << max << "]" ;
      SurfelAdjacency<2> sAdj( true );
      std::vector< std::vector< Z2i::Point >  >  vectContoursBdryPointels;
      Surfaces<Z2i::KSpace>::extractAllPointContours4C( vectContoursBdryPointels,
							ks, predicate, sAdj );  
      if(select){
	saveSelContoursAsFC(vectContoursBdryPointels,  minSize, selectCenter,  selectDistanceMax);
      }else{
	saveAllContoursAsFc(vectContoursBdryPointels,  minSize); 
      }
      trace.info() << " [done]" << endl;
    }
  }

    
  
}

