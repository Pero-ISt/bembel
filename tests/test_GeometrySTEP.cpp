#include <Bembel/Geometry>
#include "tests/Test.hpp"

bool comparePatch(Bembel::Patch &firstPatch, Bembel::Patch &secondPatch, const double epsilon = Test::Constants::test_tolerance_geometry){
    if (firstPatch.data_.size() != secondPatch.data_.size())
        {
          return false;
        }
        bool b = true;
        for (int j = 0; j < firstPatch.data_.size(); j++)
        {
          if (fabs(firstPatch.data_[j] - secondPatch.data_[j]) > epsilon) {
            std::cout << "entries are different, deviation : " << abs(firstPatch.data_[j] - secondPatch.data_[j]) << std::endl;
            std::cout <<  std::setprecision(15) << firstPatch.data_[j] << "!=" << secondPatch.data_[j] << std::endl;
            b = false;
          } 
        }
    return b;
}

/**
 * @brief runs two tests. The first tests the reading of STEP files, the second tests the writing and reading of step files
 *  
 */
int main() {
    using namespace Bembel;
    
    Geometry sphere("sphere.dat");
    
    //first Test compare STEP sphere with DAT sphere
    Geometry step_sphere("sphere.stp");

    //the STEP file has a diffrent sorting
    if(!( comparePatch(sphere.get_geometry().at(0), step_sphere.get_geometry().at(5)) &&
          comparePatch(sphere.get_geometry().at(1), step_sphere.get_geometry().at(4)) &&
          comparePatch(sphere.get_geometry().at(2), step_sphere.get_geometry().at(3)) &&
          comparePatch(sphere.get_geometry().at(3), step_sphere.get_geometry().at(2)) &&
          comparePatch(sphere.get_geometry().at(4), step_sphere.get_geometry().at(1)) &&
          comparePatch(sphere.get_geometry().at(5), step_sphere.get_geometry().at(0)))){
        std::cerr << "first Test failed" << std::endl;
    }

    //second Test export geometry form geo directory in step and import it again.
    //write step-file 
    WriteSTEPFile(sphere.get_geometry(), "sphere-step-export.step");
    //read step-file 
    Geometry geo_step("sphere-step-export.step");
        
    if(!( comparePatch(sphere.get_geometry().at(0), geo_step.get_geometry().at(0)) &&
          comparePatch(sphere.get_geometry().at(1), geo_step.get_geometry().at(1)) &&
          comparePatch(sphere.get_geometry().at(2), geo_step.get_geometry().at(2)) &&
          comparePatch(sphere.get_geometry().at(3), geo_step.get_geometry().at(3)) &&
          comparePatch(sphere.get_geometry().at(4), geo_step.get_geometry().at(4)) &&
          comparePatch(sphere.get_geometry().at(5), geo_step.get_geometry().at(5)))){
      std::cerr << "second Test failed" << std::endl;
    }
    return 0;
}
