#include "GeometryMaterialWriter.hh"

#include "G4VPhysicalVolume.hh"
#include "G4LogicalVolume.hh"
#include "G4VSolid.hh"
#include "G4Material.hh"
#include "G4SystemOfUnits.hh"
#include "G4RotationMatrix.hh"
#include "G4ThreeVector.hh"

#include <json/json.h>

#include <fstream>
#include <iostream>
#include <memory>
#include <string>

namespace
{
  Json::Value Vec3ToJsonMm(const G4ThreeVector& v)
  {
    Json::Value out(Json::arrayValue);
    out.append(v.x() / mm);
    out.append(v.y() / mm);
    out.append(v.z() / mm);
    return out;
  }

  Json::Value RotationToJson(const G4RotationMatrix& rot)
  {
    Json::Value out(Json::arrayValue);

    out.append(rot.xx()); out.append(rot.xy()); out.append(rot.xz());
    out.append(rot.yx()); out.append(rot.yy()); out.append(rot.yz());
    out.append(rot.zx()); out.append(rot.zy()); out.append(rot.zz());

    return out;
  }

  Json::Value MaterialToJson(const G4Material* mat)
  {
    Json::Value out;

    if (!mat)
    {
      out["name"] = "";
      return out;
    }

    out["name"] = mat->GetName();
    out["density_g_per_cm3"] = mat->GetDensity() / (g / cm3);
    out["radiation_length_mm"] = mat->GetRadlen() / mm;
    out["nuclear_interaction_length_mm"] = mat->GetNuclearInterLength() / mm;

    return out;
  }

  void AddPhysicalVolumeRecursive(
      const G4VPhysicalVolume* pv,
      const G4LogicalVolume* mother_lv,
      const G4ThreeVector& parent_global_translation,
      const G4RotationMatrix& parent_global_rotation,
      const std::string& parent_path,
      int depth,
      Json::Value& volumes)
  {
    if (!pv) return;

    const auto* lv = pv->GetLogicalVolume();
    auto* solid = lv ? lv->GetSolid() : nullptr;
    const auto* mat = lv ? lv->GetMaterial() : nullptr;

    const G4ThreeVector local_translation = pv->GetTranslation();

    G4RotationMatrix local_rotation;
    if (pv->GetRotation())
    {
      local_rotation = *(pv->GetRotation());
    }
    else
    {
      local_rotation = G4RotationMatrix();
    }

    const G4ThreeVector global_translation =
        parent_global_translation + parent_global_rotation * local_translation;

    const G4RotationMatrix global_rotation =
        parent_global_rotation * local_rotation;

    const std::string physical_name = pv->GetName();
    const std::string path =
        parent_path.empty() ? physical_name : parent_path + "/" + physical_name;

    Json::Value v;

    v["path"] = path;
    v["depth"] = depth;

    v["physical_name"] = physical_name;
    v["copy_no"] = pv->GetCopyNo();

    v["logical_name"] = lv ? lv->GetName() : "";
    v["mother_logical_name"] = mother_lv ? mother_lv->GetName() : "";

    v["solid_type"] = solid ? solid->GetEntityType() : "";

    if (solid)
    {
      v["cubic_volume_mm3"] = solid->GetCubicVolume() / (mm3);
      v["surface_area_mm2"] = solid->GetSurfaceArea() / (mm2);
    }
    else
    {
      v["cubic_volume_mm3"] = 0.0;
      v["surface_area_mm2"] = 0.0;
    }

    v["material"] = MaterialToJson(mat);

    v["local_translation_mm"] = Vec3ToJsonMm(local_translation);
    v["global_translation_mm"] = Vec3ToJsonMm(global_translation);

    v["local_rotation_matrix"] = RotationToJson(local_rotation);
    v["global_rotation_matrix"] = RotationToJson(global_rotation);

    volumes.append(v);

    if (!lv) return;

    const int n_daughters = lv->GetNoDaughters();

    for (int i = 0; i < n_daughters; ++i)
    {
      const auto* daughter = lv->GetDaughter(i);

      AddPhysicalVolumeRecursive(
          daughter,
          lv,
          global_translation,
          global_rotation,
          path,
          depth + 1,
          volumes);
    }
  }
}

void GeometryMaterialWriter::WriteGeometryMaterialJson(
    G4VPhysicalVolume* world,
    const std::string& output_path)
{
  Json::Value root;
  root["format"] = "cocoa_geometry_material_dump_v2";
  root["units"]["length"] = "mm";
  root["units"]["volume"] = "mm3";
  root["units"]["surface_area"] = "mm2";
  root["units"]["density"] = "g/cm3";

  root["world_physical_volume"] = world ? world->GetName() : "";

  Json::Value volumes(Json::arrayValue);

  G4ThreeVector origin(0.0, 0.0, 0.0);
  G4RotationMatrix identity;

  AddPhysicalVolumeRecursive(
      world,
      nullptr,
      origin,
      identity,
      "",
      0,
      volumes);

  root["n_physical_volumes"] = static_cast<Json::UInt64>(volumes.size());
  root["volumes"] = volumes;

  Json::StreamWriterBuilder builder;
  builder["indentation"] = "  ";

  std::ofstream out(output_path);
  if (!out)
  {
    std::cerr << "[GeometryMaterialWriter] Could not open "
              << output_path << " for writing." << std::endl;
    return;
  }

  std::unique_ptr<Json::StreamWriter> writer(builder.newStreamWriter());
  writer->write(root, &out);

  std::cout << "[GeometryMaterialWriter] Wrote geometry/material dump to "
            << output_path << std::endl;
}