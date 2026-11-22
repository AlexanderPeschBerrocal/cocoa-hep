#ifndef __GEOMETRY_MATERIAL_WRITER_HH__
#define __GEOMETRY_MATERIAL_WRITER_HH__

#include <string>

class G4Material;
class G4VPhysicalVolume;

class GeometryMaterialWriter
{
public:
  static void WriteGeometryMaterialJson(
      G4VPhysicalVolume* world,
      const std::string& output_path);

  // Mass-fraction-weighted mean atomic number and atomic mass [g/mol] of a
  // material. This is the definition ROOT uses for TGeoMixture::GetZ/GetA,
  // so the values agree with a GDML file read through TGeo.
  static void EffectiveZA(const G4Material* material, double& z_eff, double& a_eff_g_mol);
};

#endif