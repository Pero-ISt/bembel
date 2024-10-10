// This file is part of Bembel, the higher order C++ boundary element library.
//
// Copyright (C) 2024 see <http://www.bembel.eu>
//
// It was written as part of a cooperation of J. Doelz, H. Harbrecht, S. Kurz,
// M. Multerer, S. Schoeps, and F. Wolf at Technische Universitaet Darmstadt,
// Universitaet Basel, and Universita della Svizzera italiana, Lugano. This
// source code is subject to the GNU General Public License version 3 and
// provided WITHOUT ANY WARRANTY, see <http://www.bembel.eu> for further
// information.

#ifndef BEMBEL_SRC_GEOMETRY_GEOMETRYSTEP_HPP_
#define BEMBEL_SRC_GEOMETRY_GEOMETRYSTEP_HPP_


namespace Bembel {

/**
 * \ingroup Geometry
 * \brief loads geometry from file with GEOPDE-format. Note that the direction
 *        of the normals must be consistent
 *
 * \param file_name path/filename pointing to the geometry file
 * \return std::vector of NURBS::Patch describing geometry
 */
inline PatchVector LoadGeometryFileSTEP(const std::string &file_name, const Standard_CString unit = "M", const Standard_Real precision = 1.0E-15) noexcept {

  Bembel::PatchVector out;
  STEPControl_Reader reader;

  //Set Unit to METRE
  Interface_Static::SetCVal("xstep.cascade.unit", unit);
  Interface_Static::SetRVal("write.precision.val", precision);

  //read STEP file
  IFSelect_ReturnStatus status = reader.ReadFile(file_name.c_str());
  if (status != IFSelect_RetDone) {
    std::cout << "Failed to read step-file \"" << file_name << "\"." << std::endl;
  }
  reader.TransferRoots();
  TopoDS_Shape shape = reader.OneShape();

  int SurfaceCounter = 0;
  bool rational = true;
  //loop over all Faces in STEP file and extract B-spline Surface
  TopExp_Explorer explorer;
  for (explorer.Init(shape, TopAbs_FACE); explorer.More(); explorer.Next()) {
    TopoDS_Face face = TopoDS::Face(explorer.Current());
    Handle(Geom_Surface) surface = BRep_Tool::Surface(face);
    Handle(Geom_BSplineSurface) bsplineSurface = Handle(Geom_BSplineSurface)::DownCast(surface);
    if (!bsplineSurface.IsNull()) {
      SurfaceCounter++;
      Bembel::Patch tempPatch;
      std::vector<double> tempknt1;
      std::vector<double> tempknt2;
      std::vector<Eigen::MatrixXd> tmp;

      // first knotVector
      for (int i = 1; i <= bsplineSurface->NbUKnots(); i++) {
        for (int j = 0;  j < bsplineSurface->UMultiplicities()[i]; j++) {
          tempknt1.push_back(bsplineSurface->UKnot(i));
        }
      }

      // second knotVector
      for (int i = 1; i <= bsplineSurface->NbVKnots(); i++) {
        for (int j = 0;  j < bsplineSurface->VMultiplicities()[i]; j++) {
          tempknt2.push_back(bsplineSurface->VKnot(i));
        }
      }

      rational = bsplineSurface->IsURational() | bsplineSurface->IsVRational();
      //read control points and weights and write them into matrices, set weights to ones if Surface is not in NURBS format 
      //coordinate form of opencascade model: (x, y, z , w)
      //needed coordinate form for BEMBEL: (wx, wy, wz, w)
      int numUPoles = bsplineSurface->NbUPoles();
      int numVPoles = bsplineSurface->NbVPoles();
      Eigen::Matrix<double, -1, -1> tempMatrixX(numVPoles, numUPoles);
      Eigen::Matrix<double, -1, -1> tempMatrixY(numVPoles, numUPoles);
      Eigen::Matrix<double, -1, -1> tempMatrixZ(numVPoles, numUPoles);
      Eigen::Matrix<double, -1, -1> tempMatrixW(numVPoles, numUPoles);
      for (int u = 0; u < numUPoles; u++) {
        for (int v = 0; v < numVPoles; v++) {
          gp_Pnt pole = bsplineSurface->Pole(u+1, v+1);
          double weight = 1.0;
          if (rational) {
            weight = bsplineSurface->Weight(u+1, v+1); 
            tempMatrixW(v, u) = weight;  
          }
          tempMatrixX(v, u) = pole.X() * weight;
          tempMatrixY(v, u) = pole.Y() * weight;
          tempMatrixZ(v, u) = pole.Z() * weight;
        }   
      }
      if (!rational)
      {
        tempMatrixW.setOnes();
      }

      //Init matricies in Bembel Patch
      tmp.push_back(tempMatrixX);
      tmp.push_back(tempMatrixY);
      tmp.push_back(tempMatrixZ);
      tmp.push_back(tempMatrixW);
      tempPatch.init_Patch(tmp, tempknt1, tempknt2);
      out.push_back(tempPatch);
    }
  }
   auto bspline = rational ? " NURBS":" B-spline";
  std::cout << SurfaceCounter << bspline <<" surfaces were extracted from file \"" << file_name << "\"." << std::endl;
  return out;
}

/**
 * \ingroup Geometry
 * \brief exports a geometry from Bembel to a .stp file.
 *
 * This functions assumes a p-open knot vector without internal
 * knots.
 * 
 * \param Geometry std::vector of NURBS::Patch describing geometry
 * \param file_name path/filename to be written to
 */
void WriteSTEPFile(const std::vector<Patch> &Geometry, const std::string &file_name, const Standard_CString unit = "M", const Standard_Real precision = 1.0E-15) {


  //Define opencascade compound object
  TopoDS_Compound compound;
  BRep_Builder builder;
  builder.MakeCompound(compound);

  for (auto &patch : Geometry){
    //Controlpoint and weight Arrays
    TColgp_Array2OfPnt controlPoints(1, patch.polynomial_degree_x_, 1, patch.polynomial_degree_y_);
    TColStd_Array2OfReal weights(1, patch.polynomial_degree_x_, 1, patch.polynomial_degree_y_);

    // Degree of NURBS-Surface in U- and V-direction
    Standard_Integer degreeU = patch.polynomial_degree_x_ - 1;
    Standard_Integer degreeV = patch.polynomial_degree_y_ - 1;

    // Knots and Mults
    TColStd_Array1OfReal knotsU(1, 2);
    knotsU.SetValue(1, 0.0);
    knotsU.SetValue(2, 1.0);

    TColStd_Array1OfInteger multsU(1, 2);
    multsU.SetValue(1, patch.polynomial_degree_x_);
    multsU.SetValue(2, patch.polynomial_degree_x_);

    TColStd_Array1OfReal knotsV(1, 2);
    knotsV.SetValue(1, 0.0);
    knotsV.SetValue(2, 1.0);

    TColStd_Array1OfInteger multsV(1, 2);
    multsV.SetValue(1, patch.polynomial_degree_y_);
    multsV.SetValue(2, patch.polynomial_degree_y_);

    //set control-points and weights in form x,y,z,w (Bembel form: wx, wy, wz, w)
    int index = 0;
    for (int i = 1; i < patch.polynomial_degree_x_+1; i++)
    {
      for (int j = 1; j < patch.polynomial_degree_y_+1; j++)
      {
      double weight = patch.data_[index+3];
      double x = patch.data_[index] / weight;
      double y = patch.data_[index+1] / weight;
      double z = patch.data_[index+2] / weight;
      controlPoints.SetValue(i, j, gp_Pnt(x, y, z));
      weights.SetValue(i, j, weight);
      index+=4;
      }
    }

    //Define NURBS-Surface
    Handle(Geom_BSplineSurface) nurbsSurface = new Geom_BSplineSurface(
        controlPoints,
        weights,
        knotsU,
        knotsV,
        multsU,
        multsV,
        degreeU,
        degreeV
    );
    //Add Surface to Compound Object
    builder.Add(compound, BRepBuilderAPI_MakeFace(nurbsSurface, 0).Face());
  }

  //Set Unit to METRE + Set precision
  Interface_Static::SetCVal("xstep.cascade.unit", unit);
  Interface_Static::SetRVal("write.precision.val", precision);
  
  //Write Step file
  STEPControl_Writer writer;
  writer.Transfer(compound, STEPControl_AsIs);
  IFSelect_ReturnStatus status = writer.Write(file_name.c_str());

  if (status != IFSelect_RetDone) {
    std::cerr << "Failed to write step-file \"" << file_name << "\"." << std::endl;
  }

  return;
}

} // namespace Bembel
#endif  // BEMBEL_SRC_GEOMETRY_GEOMETRYSTEP_HPP_
