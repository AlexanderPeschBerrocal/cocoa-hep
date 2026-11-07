#ifndef __GEOMETRY_MATERIAL_WRITER_HH__
#define __GEOMETRY_MATERIAL_WRITER_HH__

#include <string>

class G4VPhysicalVolume;

class GeometryMaterialWriter
{
public:
  static void WriteGeometryMaterialJson(
      G4VPhysicalVolume* world,
      const std::string& output_path);
};

#endif