#include "InnerConstruction.hh"

#include "G4Box.hh"
#include "G4ChordFinder.hh"
#include "G4Element.hh"
#include "G4NistManager.hh"
#include "G4FieldManager.hh"
#include "G4LogicalVolume.hh"

#include "G4PVPlacement.hh"
#include "G4SDManager.hh"
#include "G4SystemOfUnits.hh"
#include "G4TransportationManager.hh"
#include "G4Tubs.hh"
#include "G4CutTubs.hh"
#include "G4Cons.hh"

#include "G4SystemOfUnits.hh"
#include "G4IntersectionSolid.hh"
#include "G4PVPlacement.hh"
#include "G4SubtractionSolid.hh"
#include "G4PVParameterised.hh"
#include "DetectorGeometryDefinitions.hh"
#include "G4ReflectionFactory.hh"
#include "G4UniformMagField.hh"
#include "G4FieldManager.hh"

#include "G4RunManager.hh"

#include <array>
#include <cmath>
#include <string>

#include "G4Region.hh"

InnerConstruction::InnerConstruction(G4LogicalVolume* expHallLV, G4Material* default_Material, G4Material* Iron, G4Material* ElSi, bool fCheck_Overlaps) 
{

    //
    // All inner detector logical volume names start with "Inner".
    // This property will be used elsewhere.
    //

    theta_min = 2*atan(exp(-1*config_obj.max_eta_barrel));
    fCheckOverlaps = fCheck_Overlaps;
    Iron_Support_VisAtt->SetForceSolid(true);
    defaultMaterial = default_Material;
    r_inn = config_obj.r_inn_calo;
    iron = Iron;
    elSi = ElSi;



    //* Inner created to have magnetic field only in the inner part of detector
	G4Tubs *Inner= new G4Tubs("Inner_det", 0, r_inn, r_inn/tan(theta_min), 0, config_obj.max_phi);
	GlobalLV = new G4LogicalVolume(Inner, defaultMaterial, "Inner_LV");
	new G4PVPlacement(
					0,                // no rotation
					G4ThreeVector(0., 0., 0. ),
					GlobalLV,          // its logical volume
					"Inner_PL",          // its name
					expHallLV,                // its mother  volume
					false,            // no boolean operation
					0  ,              // copy number
					fCheckOverlaps
	);

	float       SpaceForLayers        = r_inn - ( r_out_trkStr3 + widthiron_add );
	float       GapBetweenIronLayers  = SpaceForLayers/(NumberOfIronLayers+1);
	long double r_magField            = (r_out_trkStr3 + widthiron_add) + GapBetweenIronLayers;
	long double iron_r_inn_ec         = ( r_out_trkStr3 + widthiron_add ) + NumberOfIronLayers * GapBetweenIronLayers;
	long double iron_width_ec         = 4.4 * cm;
	long double l_magField            = ( ( iron_r_inn_ec +0.5 * iron_width_ec ) / tan( theta_min ) ); // gap offset
	l_magField                       -= 0.5 * fabs( iron_width_ec / ( tan( theta_min ) ) );            // gap half-width

	G4Tubs *magFieldTub = new G4Tubs( "Magnetic Field Tube",
					  0.0,
					  r_magField,
					  l_magField,
					  0.0,
					  config_obj.max_phi );
	MagFieldLV = new G4LogicalVolume( magFieldTub,
					  defaultMaterial,
					  "MagField_LV" );
	new G4PVPlacement( 0,
			   G4ThreeVector(0., 0., 0. ),
			   MagFieldLV,
			   "Inner_MagField_PL",
			   GlobalLV,
			   false,
			   0,
			   fCheckOverlaps );
	
	  Barrel_Inner();
	  EndCap_Inner();
	  if ( config_obj.use_inner_detector && config_obj.use_detailed_tracker_material ) {
	      BuildDetailedTrackerMaterial();
	  }
	
	auto* trackerRegion = new G4Region("TrackerRegion");
	GlobalLV->SetRegion(trackerRegion);
	trackerRegion->AddRootLogicalVolume(GlobalLV);

	if (config_obj.fieldValue!=0)
	{
		G4UniformMagField* magField = new G4UniformMagField(G4ThreeVector(0.,0.,config_obj.fieldValue));
		G4FieldManager* fieldMgr= new G4FieldManager(magField);
		MagFieldLV->SetFieldManager(fieldMgr,true);
		fieldMgr->SetDetectorField(magField);
		fieldMgr->CreateChordFinder(magField);
		fieldMgr->GetChordFinder()->SetDeltaChord(0.01);
	}

}
InnerConstruction::~InnerConstruction()
{;}

void InnerConstruction::BuildDetailedTrackerMaterial()
{
    // This is a deliberately simplified passive-material model.  The dimensions
    // leave clearance around the idealised silicon and iron tracker layers.
    auto* nist = G4NistManager::Instance();
    auto* carbon = nist->FindOrBuildMaterial("G4_C");
    auto* aluminium = nist->FindOrBuildMaterial("G4_Al");
    auto* coolant = nist->FindOrBuildMaterial("G4_WATER");
    auto* circuitBoard = nist->FindOrBuildMaterial("G4_POLYETHYLENE");
    auto* copper = nist->FindOrBuildMaterial("G4_Cu");
    auto* cableMaterial = nist->FindOrBuildMaterial("G4_KAPTON");

    auto* carbonVis = new G4VisAttributes(true, G4Colour(0.25, 0.25, 0.25));
    auto* pipeVis = new G4VisAttributes(true, G4Colour(0.65, 0.65, 0.65));
    auto* coolantVis = new G4VisAttributes(true, G4Colour(0.1, 0.7, 1.0));
    auto* electronicsVis = new G4VisAttributes(true, G4Colour(0.1, 0.55, 0.1));
    auto* copperVis = new G4VisAttributes(true, G4Colour(0.8, 0.4, 0.1));
    auto* cableVis = new G4VisAttributes(true, G4Colour(0.8, 0.65, 0.1));
    carbonVis->SetForceSolid(true);
    pipeVis->SetForceSolid(true);
    coolantVis->SetForceSolid(true);
    electronicsVis->SetForceSolid(true);
    copperVis->SetForceSolid(true);
    cableVis->SetForceSolid(true);

    struct BarrelLayer {
        G4double outerRadius;
        G4double halfLength;
    };
    const std::array<BarrelLayer, 9> barrelLayers = {{
        {r_out_trkPix0, 280 * mm}, {r_out_trkPix1, 280 * mm},
        {r_out_trkPix2, 280 * mm}, {r_out_trkPix3, 280 * mm},
        {r_out_trkPix4, 280 * mm}, {r_out_trkStr0, 1150 * mm},
        {r_out_trkStr1, 1150 * mm}, {r_out_trkStr2, 1150 * mm},
        {r_out_trkStr3, 1150 * mm}
    }};

    constexpr G4int pipesPerLayer = 8;
    const G4double carbonThickness = 0.5 * mm;
    const G4double pipeOuterRadius = 1.5 * mm;
    const G4double pipeInnerRadius = 1.1 * mm;

    for (std::size_t layerIndex = 0; layerIndex < barrelLayers.size(); ++layerIndex) {
        const auto& layer = barrelLayers[layerIndex];
        const std::string suffix = std::to_string(layerIndex);
        const G4double carbonInnerRadius = layer.outerRadius + widthiron_add + 0.2 * mm;
        const G4double serviceRadius = layer.outerRadius + widthiron_add + 3.0 * mm;

        auto* supportSolid = new G4Tubs(
            "Inner_Detailed_Carbon_Solid_" + suffix,
            carbonInnerRadius,
            carbonInnerRadius + carbonThickness,
            layer.halfLength,
            0.,
            config_obj.max_phi);
        auto* supportLV = new G4LogicalVolume(
            supportSolid, carbon, "Inner_Detailed_Carbon_LV_" + suffix);
        new G4PVPlacement(
            nullptr, G4ThreeVector(), supportLV,
            "Inner_Detailed_Carbon_PL_" + suffix,
            MagFieldLV, false, static_cast<G4int>(layerIndex), fCheckOverlaps);
        supportLV->SetVisAttributes(carbonVis);

        auto* pipeSolid = new G4Tubs(
            "Inner_Detailed_CoolingPipe_Solid_" + suffix,
            0., pipeOuterRadius, layer.halfLength, 0., config_obj.max_phi);
        auto* pipeLV = new G4LogicalVolume(
            pipeSolid, aluminium, "Inner_Detailed_CoolingPipe_LV_" + suffix);
        auto* coolantSolid = new G4Tubs(
            "Inner_Detailed_Coolant_Solid_" + suffix,
            0., pipeInnerRadius, layer.halfLength - 0.1 * mm,
            0., config_obj.max_phi);
        auto* coolantLV = new G4LogicalVolume(
            coolantSolid, coolant, "Inner_Detailed_Coolant_LV_" + suffix);
        new G4PVPlacement(
            nullptr, G4ThreeVector(), coolantLV,
            "Inner_Detailed_Coolant_PL_" + suffix,
            pipeLV, false, 0, fCheckOverlaps);
        pipeLV->SetVisAttributes(pipeVis);
        coolantLV->SetVisAttributes(coolantVis);

        for (G4int pipeIndex = 0; pipeIndex < pipesPerLayer; ++pipeIndex) {
            const G4double phi = config_obj.max_phi * pipeIndex / pipesPerLayer;
            new G4PVPlacement(
                nullptr,
                G4ThreeVector(serviceRadius * std::cos(phi),
                              serviceRadius * std::sin(phi), 0.),
                pipeLV,
                "Inner_Detailed_CoolingPipe_PL_" + suffix + "_" + std::to_string(pipeIndex),
                MagFieldLV, false, pipeIndex, fCheckOverlaps);
        }

        // Polymer-substrate and copper annuli model front-end boards at both barrel ends.
        const G4double boardInnerRadius = layer.outerRadius - 1.0 * mm;
        const G4double boardOuterRadius = serviceRadius + 2.0 * mm;
        auto* boardSolid = new G4Tubs(
            "Inner_Detailed_ElectronicsBoard_Solid_" + suffix,
            boardInnerRadius, boardOuterRadius, 1.0 * mm, 0., config_obj.max_phi);
        auto* boardLV = new G4LogicalVolume(
            boardSolid, circuitBoard, "Inner_Detailed_ElectronicsBoard_LV_" + suffix);
        auto* copperSolid = new G4Tubs(
            "Inner_Detailed_ElectronicsCopper_Solid_" + suffix,
            boardInnerRadius, boardOuterRadius, 0.1 * mm, 0., config_obj.max_phi);
        auto* copperLV = new G4LogicalVolume(
            copperSolid, copper, "Inner_Detailed_ElectronicsCopper_LV_" + suffix);
        boardLV->SetVisAttributes(electronicsVis);
        copperLV->SetVisAttributes(copperVis);

        for (G4int direction : {-1, 1}) {
            new G4PVPlacement(
                nullptr, G4ThreeVector(0., 0., direction * (layer.halfLength + 3.0 * mm)),
                boardLV,
                "Inner_Detailed_ElectronicsBoard_PL_" + suffix + "_" + std::to_string(direction),
                MagFieldLV, false, direction, fCheckOverlaps);
            new G4PVPlacement(
                nullptr, G4ThreeVector(0., 0., direction * (layer.halfLength + 4.2 * mm)),
                copperLV,
                "Inner_Detailed_ElectronicsCopper_PL_" + suffix + "_" + std::to_string(direction),
                MagFieldLV, false, direction, fCheckOverlaps);
        }
    }

    // Carbon-fibre backing disks follow each silicon endcap disk.
    struct EndcapLayer {
        G4double innerRadius;
        G4double outerRadius;
        G4double position;
    };
    const std::array<EndcapLayer, 15> endcapLayers = {{
        {r_inn_trkPix0, r_out_trkPix4, pos_EndCap_trkPix0},
        {r_inn_trkPix0, r_out_trkPix4, pos_EndCap_trkPix1},
        {r_inn_trkPix0, r_out_trkPix4, pos_EndCap_trkPix2},
        {r_inn_trkPix0, r_out_trkPix4, pos_EndCap_trkPix3},
        {r_inn_trkPix0, r_out_trkPix4, pos_EndCap_trkPix4},
        {r_inn_trkPix0, r_out_trkPix4, pos_EndCap_trkPix5},
        {r_inn_trkPix0, r_out_trkPix4, pos_EndCap_trkPix6},
        {r_inn_trkPix2, r_out_trkPix4, pos_EndCap_trkPix7},
        {r_inn_trkPix2, r_out_trkPix4, pos_EndCap_trkPix8},
        {r_inn_trkPix2, r_out_trkPix4, pos_EndCap_trkPix9},
        {r_inn_trkStr0, r_out_trkStr3, pos_EndCap_trkStr0},
        {r_inn_trkStr0, r_out_trkStr3, pos_EndCap_trkStr1},
        {r_inn_trkStr0, r_out_trkStr3, pos_EndCap_trkStr2},
        {r_inn_trkStr0, r_out_trkStr3, pos_EndCap_trkStr3},
        {r_inn_trkStr0, r_out_trkStr3, pos_EndCap_trkStr4}
    }};

    for (std::size_t layerIndex = 0; layerIndex < endcapLayers.size(); ++layerIndex) {
        const auto& layer = endcapLayers[layerIndex];
        const std::string suffix = std::to_string(layerIndex);
        auto* backingSolid = new G4Tubs(
            "Inner_Detailed_EndcapBacking_Solid_" + suffix,
            layer.innerRadius, layer.outerRadius, 0.25 * mm, 0., config_obj.max_phi);
        auto* backingLV = new G4LogicalVolume(
            backingSolid, carbon, "Inner_Detailed_EndcapBacking_LV_" + suffix);
        backingLV->SetVisAttributes(carbonVis);
        for (G4int direction : {-1, 1}) {
            new G4PVPlacement(
                nullptr,
                G4ThreeVector(0., 0., direction * (layer.position + widthiron_add + 0.8 * mm)),
                backingLV,
                "Inner_Detailed_EndcapBacking_PL_" + suffix + "_" + std::to_string(direction),
                MagFieldLV, false, direction, fCheckOverlaps);
        }
    }

    // A thin service/cable cylinder joins the strip barrel to two endplates.
    auto* cableSolid = new G4Tubs(
        "Inner_Detailed_CableBundle_Solid", 1010 * mm, 1012 * mm,
        1250 * mm, 0., config_obj.max_phi);
    auto* cableLV = new G4LogicalVolume(
        cableSolid, cableMaterial, "Inner_Detailed_CableBundle_LV");
    new G4PVPlacement(
        nullptr, G4ThreeVector(), cableLV, "Inner_Detailed_CableBundle_PL",
        MagFieldLV, false, 0, fCheckOverlaps);
    cableLV->SetVisAttributes(cableVis);

    auto* endplateSolid = new G4Tubs(
        "Inner_Detailed_ServiceEndplate_Solid", 30 * mm, 1050 * mm,
        2.0 * mm, 0., config_obj.max_phi);
    auto* endplateLV = new G4LogicalVolume(
        endplateSolid, circuitBoard, "Inner_Detailed_ServiceEndplate_LV");
    auto* endplateCopperSolid = new G4Tubs(
        "Inner_Detailed_ServiceEndplateCopper_Solid", 30 * mm, 1050 * mm,
        0.1 * mm, 0., config_obj.max_phi);
    auto* endplateCopperLV = new G4LogicalVolume(
        endplateCopperSolid, copper, "Inner_Detailed_ServiceEndplateCopper_LV");
    endplateLV->SetVisAttributes(electronicsVis);
    endplateCopperLV->SetVisAttributes(copperVis);
    for (G4int direction : {-1, 1}) {
        new G4PVPlacement(
            nullptr, G4ThreeVector(0., 0., direction * 2700 * mm),
            endplateLV,
            "Inner_Detailed_ServiceEndplate_PL_" + std::to_string(direction),
            MagFieldLV, false, direction, fCheckOverlaps);
        new G4PVPlacement(
            nullptr, G4ThreeVector(0., 0., direction * 2702.25 * mm),
            endplateCopperLV,
            "Inner_Detailed_ServiceEndplateCopper_PL_" + std::to_string(direction),
            MagFieldLV, false, direction, fCheckOverlaps);
    }
}

void InnerConstruction::Barrel_Inner()
{
	if ( config_obj.use_inner_detector ) {
	
	    PixelTrk_Barrel(r_inn_trkPix0,r_out_trkPix0,280*mm,Pix_VisAtt);
	    PixelTrk_Barrel(r_inn_trkPix1,r_out_trkPix1,280*mm,Pix_VisAtt);
	    PixelTrk_Barrel(r_inn_trkPix2,r_out_trkPix2,280*mm,Pix_VisAtt);
	    PixelTrk_Barrel(r_inn_trkPix3,r_out_trkPix3,280*mm,Pix_VisAtt);
	    PixelTrk_Barrel(r_inn_trkPix4,r_out_trkPix4,280*mm,Pix_VisAtt);
	    
	    
	    PixelTrk_Barrel(r_inn_trkStr0,r_out_trkStr0,1150*mm,Str_VisAtt);
	    PixelTrk_Barrel(r_inn_trkStr1,r_out_trkStr1,1150*mm,Str_VisAtt);
	    PixelTrk_Barrel(r_inn_trkStr2,r_out_trkStr2,1150*mm,Str_VisAtt);
	    PixelTrk_Barrel(r_inn_trkStr3,r_out_trkStr3,1150*mm,Str_VisAtt, true);

	}

	if ( !config_obj.use_ID_support )
	    return;
	    
    // *Barrel support material 
	float SpaceForLayers = r_inn-(r_out_trkStr3+widthiron_add);
	float GapBetweenIronLayers = SpaceForLayers/(NumberOfIronLayers+1);
	long double iron_width_barrel = 4.4*cm/NumberOfIronLayers;
	G4LogicalVolume *Iron_Layer_LV;
	G4Tubs *Iron_Layer;
	for (int NironLayer = 1; NironLayer < (NumberOfIronLayers+1); NironLayer++)
	{

		long double iron_r_inn = (r_out_trkStr3 + widthiron_add) + NironLayer*GapBetweenIronLayers;
		
		long double l_Iron_Layer =  (iron_r_inn)/tan(theta_min);

		if (NironLayer!=NumberOfIronLayers)
		{
			l_Iron_Layer =  (iron_r_inn+0.9*iron_width_barrel)/tan(theta_min);
			Iron_Layer = new G4Tubs("Inner_Iron_layer", iron_r_inn, iron_r_inn + iron_width_barrel, l_Iron_Layer, 0, config_obj.max_phi);
			Iron_Layer_LV = new G4LogicalVolume(Iron_Layer, iron, "Iron_gap_barrel_ID_LV");
			new G4PVPlacement(
							0,                // no rotation
							G4ThreeVector(0., 0., 0. ),
							Iron_Layer_LV,          // its logical volume
							( std::string( "Iron_gap_barrel_" ) + std::to_string( NironLayer ) ).c_str(),          // its name
							GlobalLV,                // its mother  volume
							false,            // no boolean operation
							0,               // copy number
							fCheckOverlaps
						);
		}
		else 
		{
			Iron_Layer = new G4Tubs("Inner_Iron_layer", iron_r_inn, iron_r_inn + iron_width_barrel, l_Iron_Layer, 0, config_obj.max_phi);
			Iron_Layer_LV = new G4LogicalVolume(Iron_Layer, iron, "Iron_gap_barrel_ID_LV");
			new G4PVPlacement(
								0,                // no rotation
								G4ThreeVector(0., 0., 0. ),
								Iron_Layer_LV,          // its logical volume
								( std::string( "Iron_gap_barrel_" ) + std::to_string( NironLayer ) ).c_str(),          // its name
								GlobalLV,                // its mother  volume
								false,            // no boolean operation
								0,               // copy number
								fCheckOverlaps
							);
		}
		
		Iron_Layer_LV->SetVisAttributes(Iron_Support_VisAtt);
		
	}
}

void InnerConstruction::EndCap_Inner()
{

    if ( config_obj.use_inner_detector ) {
	
	PixelTrk_EndCap(r_inn_trkPix0,r_out_trkPix4,r_out_trkPix0-r_inn_trkPix0,pos_EndCap_trkPix0,Pix_VisAtt, 1);
	PixelTrk_EndCap(r_inn_trkPix0,r_out_trkPix4,r_out_trkPix0-r_inn_trkPix0,pos_EndCap_trkPix0,Pix_VisAtt,-1);
	
	PixelTrk_EndCap(r_inn_trkPix0,r_out_trkPix4,r_out_trkPix0-r_inn_trkPix0,pos_EndCap_trkPix1,Pix_VisAtt, 1);
	PixelTrk_EndCap(r_inn_trkPix0,r_out_trkPix4,r_out_trkPix0-r_inn_trkPix0,pos_EndCap_trkPix1,Pix_VisAtt,-1);
	
	PixelTrk_EndCap(r_inn_trkPix0,r_out_trkPix4,r_out_trkPix0-r_inn_trkPix0,pos_EndCap_trkPix2,Pix_VisAtt, 1);
	PixelTrk_EndCap(r_inn_trkPix0,r_out_trkPix4,r_out_trkPix0-r_inn_trkPix0,pos_EndCap_trkPix2,Pix_VisAtt,-1);
	
	PixelTrk_EndCap(r_inn_trkPix0,r_out_trkPix4,r_out_trkPix0-r_inn_trkPix0,pos_EndCap_trkPix3,Pix_VisAtt, 1);
	PixelTrk_EndCap(r_inn_trkPix0,r_out_trkPix4,r_out_trkPix0-r_inn_trkPix0,pos_EndCap_trkPix3,Pix_VisAtt,-1);
	
	PixelTrk_EndCap(r_inn_trkPix0,r_out_trkPix4,r_out_trkPix0-r_inn_trkPix0,pos_EndCap_trkPix4,Pix_VisAtt, 1);
	PixelTrk_EndCap(r_inn_trkPix0,r_out_trkPix4,r_out_trkPix0-r_inn_trkPix0,pos_EndCap_trkPix4,Pix_VisAtt,-1);
	
	PixelTrk_EndCap(r_inn_trkPix0,r_out_trkPix4,r_out_trkPix0-r_inn_trkPix0,pos_EndCap_trkPix5,Pix_VisAtt, 1);
	PixelTrk_EndCap(r_inn_trkPix0,r_out_trkPix4,r_out_trkPix0-r_inn_trkPix0,pos_EndCap_trkPix5,Pix_VisAtt,-1);
	
	PixelTrk_EndCap(r_inn_trkPix0,r_out_trkPix4,r_out_trkPix0-r_inn_trkPix0,pos_EndCap_trkPix6,Pix_VisAtt, 1);
	PixelTrk_EndCap(r_inn_trkPix0,r_out_trkPix4,r_out_trkPix0-r_inn_trkPix0,pos_EndCap_trkPix6,Pix_VisAtt,-1);

	PixelTrk_EndCap(r_inn_trkPix2,r_out_trkPix4,r_out_trkPix0-r_inn_trkPix0,pos_EndCap_trkPix7,Pix_VisAtt, 1);
	PixelTrk_EndCap(r_inn_trkPix2,r_out_trkPix4,r_out_trkPix0-r_inn_trkPix0,pos_EndCap_trkPix7,Pix_VisAtt,-1);
	
	PixelTrk_EndCap(r_inn_trkPix2,r_out_trkPix4,r_out_trkPix0-r_inn_trkPix0,pos_EndCap_trkPix8,Pix_VisAtt, 1);
	PixelTrk_EndCap(r_inn_trkPix2,r_out_trkPix4,r_out_trkPix0-r_inn_trkPix0,pos_EndCap_trkPix8,Pix_VisAtt,-1);

	PixelTrk_EndCap(r_inn_trkPix2,r_out_trkPix4,r_out_trkPix0-r_inn_trkPix0,pos_EndCap_trkPix9,Pix_VisAtt, 1);
	PixelTrk_EndCap(r_inn_trkPix2,r_out_trkPix4,r_out_trkPix0-r_inn_trkPix0,pos_EndCap_trkPix9,Pix_VisAtt,-1);
	
	
	PixelTrk_EndCap(r_inn_trkStr0,r_out_trkStr3,r_out_trkPix0-r_inn_trkPix0,pos_EndCap_trkStr0,Str_VisAtt, 1);
	PixelTrk_EndCap(r_inn_trkStr0,r_out_trkStr3,r_out_trkPix0-r_inn_trkPix0,pos_EndCap_trkStr0,Str_VisAtt,-1);
	
	PixelTrk_EndCap(r_inn_trkStr0,r_out_trkStr3,r_out_trkPix0-r_inn_trkPix0,pos_EndCap_trkStr1,Str_VisAtt, 1);
	PixelTrk_EndCap(r_inn_trkStr0,r_out_trkStr3,r_out_trkPix0-r_inn_trkPix0,pos_EndCap_trkStr1,Str_VisAtt,-1);
	
	PixelTrk_EndCap(r_inn_trkStr0,r_out_trkStr3,r_out_trkPix0-r_inn_trkPix0,pos_EndCap_trkStr2,Str_VisAtt, 1);
	PixelTrk_EndCap(r_inn_trkStr0,r_out_trkStr3,r_out_trkPix0-r_inn_trkPix0,pos_EndCap_trkStr2,Str_VisAtt,-1);

	PixelTrk_EndCap(r_inn_trkStr0,r_out_trkStr3,r_out_trkPix0-r_inn_trkPix0,pos_EndCap_trkStr3,Str_VisAtt, 1);
	PixelTrk_EndCap(r_inn_trkStr0,r_out_trkStr3,r_out_trkPix0-r_inn_trkPix0,pos_EndCap_trkStr3,Str_VisAtt,-1);

	PixelTrk_EndCap(r_inn_trkStr0,r_out_trkStr3,r_out_trkPix0-r_inn_trkPix0,pos_EndCap_trkStr4,Str_VisAtt, 1, true);
	PixelTrk_EndCap(r_inn_trkStr0,r_out_trkStr3,r_out_trkPix0-r_inn_trkPix0,pos_EndCap_trkStr4,Str_VisAtt,-1, true);

    }
    
	if ( !config_obj.use_ID_support ) return;

    
	//* End-Cap support material 
	float SpaceForLayers = r_inn-(r_out_trkStr3+widthiron_add);
	float GapBetweenIronLayers = SpaceForLayers/(NumberOfIronLayers+1);
	long double iron_r_inn = (r_out_trkStr3 + widthiron_add) + NumberOfIronLayers*GapBetweenIronLayers;
	long double iron_width_endcap = 4.4*cm;
	long double iron_width_barrel = 4.4*cm/NumberOfIronLayers;

	for (int direction = 1; direction > -2; direction=direction-2)
	{
		long double length_cone_min = iron_r_inn/tan(theta_min);
		long double length_cone_max = length_cone_min+iron_width_endcap/tan(theta_min);
		long double d_theta_next = 2*atan(exp(-1*(config_obj.max_eta_endcap)));
		if (direction==-1)
		{
			long double buf = length_cone_min;
			length_cone_min = length_cone_max;
			length_cone_max = buf;
		}

		G4Cons *Cone_Trk = new G4Cons("Cone_Gap", (length_cone_max)*tan(d_theta_next), (iron_r_inn + iron_width_barrel) , (length_cone_min)*tan(d_theta_next), (iron_r_inn + iron_width_barrel), fabs(iron_width_endcap/(2*tan(theta_min))),  0, 2*M_PI);

		G4LogicalVolume *Iron_LV_posdir = new G4LogicalVolume(Cone_Trk, iron, "Iron_ID_cone_endcap");

		new G4PVPlacement(
						0,               //* rotation
						G4ThreeVector(0, 0,  -1* direction*((iron_r_inn+0.5*iron_width_endcap)/tan(theta_min))), //* Placed the pixel in specific way, that its slices goes throught (0 ,0, 0)
						Iron_LV_posdir,          //* its logical volume
						"Iron_ID_gap_endcap",          //* its name
						GlobalLV,                //* its mother  volume
						false,            //* no boolean operation
						0,             //* copy number
						fCheckOverlaps
							);
		Iron_LV_posdir->SetVisAttributes(Iron_Support_VisAtt);
	}



}


void InnerConstruction::PixelTrk_EndCap(long double r_inn_trkPix, long double r_out_trkPix, G4VisAttributes* VisAtt, int direction) 
{
	
	
	long double length_cone_min = (r_inn_trkPix)/tan(theta_min);
	long double length_cone_max = length_cone_min+(r_out_trkPix-r_inn_trkPix)/tan(theta_min);
	if (direction == -1)
	{
		long double buf = length_cone_min;
		length_cone_min = length_cone_max;
		length_cone_max = buf;
	}
	long double d_theta = 2*atan(exp(-1*config_obj.max_eta_barrel));
	long double d_theta_next = 2*atan(exp(-1*(config_obj.max_eta_endcap)));


	G4Cons *Cone_Trk = new G4Cons("Cone_Gap", (length_cone_max)*tan(d_theta_next), length_cone_max*tan(d_theta),(length_cone_min)*tan(d_theta_next),length_cone_min*tan(d_theta), fabs((length_cone_max-length_cone_min)/2),  0, 2*M_PI);

	G4LogicalVolume *Trk_LV_posdir = new G4LogicalVolume(Cone_Trk, elSi, "Inner_Cone_Gap");

	new G4PVPlacement(
					0,               //* rotation
					G4ThreeVector(0, 0,  -1* direction*(length_cone_max+length_cone_min)/2), //* Placed the pixel in specific way, that its slices goes throught (0 ,0, 0)
					Trk_LV_posdir,          //* its logical volume
					"Iron_PL",          //* its name
					MagFieldLV,                //* its mother  volume
					false,            //* no boolean operation
					0,             //* copy number
					fCheckOverlaps
		);


	length_cone_min = (r_out_trkPix)/tan(theta_min);
	length_cone_max = length_cone_min+(widthiron_add)/tan(theta_min);
	if (direction == -1)
	{
		long double buf = length_cone_min;
		length_cone_min = length_cone_max;
		length_cone_max = buf;
	}
	d_theta = 2*atan(exp(-1*config_obj.max_eta_barrel));
	d_theta_next = 2*atan(exp(-1*config_obj.max_eta_endcap));

	G4Cons *Cone_Trk_iron = new G4Cons("Cone_Gap", (length_cone_max)*tan(d_theta_next), length_cone_max*tan(d_theta),(length_cone_min)*tan(d_theta_next),length_cone_min*tan(d_theta), fabs((length_cone_max-length_cone_min)/2),  0, 2*M_PI);
	G4LogicalVolume *Inner_trkPix_LV_iron = new G4LogicalVolume(Cone_Trk_iron, iron, "Inner_trkPix_LV_iron");
	new G4PVPlacement(
					0,               //* rotation
					G4ThreeVector(0, 0,  -1* direction*(length_cone_max+length_cone_min)/2), //* Placed the pixel in specific way, that its slices goes throught (0 ,0, 0)
					Inner_trkPix_LV_iron,          //* its logical volume
					"Iron_PL",          //* its name
					MagFieldLV,                //* its mother  volume
					false,            //* no boolean operation
					0,             //* copy number
					fCheckOverlaps
		);

	Trk_LV_posdir->SetVisAttributes(VisAtt);
	VisAtt->SetForceSolid(true);
	Inner_trkPix_LV_iron->SetVisAttributes(VisAtt);
}
void InnerConstruction::PixelTrk_EndCap(long double r_inn_trkPix, long double r_out_trkPix, long double width, long double Abs_Position, G4VisAttributes* VisAtt, int direction, bool isOutermostLayer) 
{


	G4Tubs *Cone_Trk = new G4Tubs("Cone_Gap", r_inn_trkPix, r_out_trkPix, width/2,  0, 2*M_PI);

	G4LogicalVolume *Trk_LV_posdir = new G4LogicalVolume(Cone_Trk, elSi, "Inner_Cone_Gap");

	std::string detName = "Inner_si_endcap_PL";
	std::string supName = "Inner_Iron_PL";
	if ( isOutermostLayer ) {
	    detName += "_outermostInner";
	    supName += "_outermostInner";
	}
	
	new G4PVPlacement(
					0,               //* rotation
					G4ThreeVector(0, 0,  -1* direction*Abs_Position), //* Placed the pixel in specific way, that its slices goes throught (0 ,0, 0)
					Trk_LV_posdir,          //* its logical volume
					detName,          //* its name
					MagFieldLV,                //* its mother  volume
					false,            //* no boolean operation
					0,             //* copy number
					fCheckOverlaps
		);


	G4Tubs *Cone_Trk_iron = new G4Tubs("Cone_Gap", r_inn_trkPix, r_out_trkPix, widthiron_add/2,  0, 2*M_PI);
	G4LogicalVolume *Inner_trkPix_LV_iron = new G4LogicalVolume(Cone_Trk_iron, iron, "Inner_trkPix_LV_iron");
	new G4PVPlacement(
					0,               //* rotation
					G4ThreeVector(0, 0,  -1* direction*(Abs_Position+widthiron_add)), //* Placed the pixel in specific way, that its slices goes throught (0 ,0, 0)
					Inner_trkPix_LV_iron,          //* its logical volume
					supName,          //* its name
					MagFieldLV,                //* its mother  volume
					false,            //* no boolean operation
					0,             //* copy number
					fCheckOverlaps
		);

	Trk_LV_posdir->SetVisAttributes(VisAtt);
	VisAtt->SetForceSolid(true);
	Inner_trkPix_LV_iron->SetVisAttributes(VisAtt);
}

void InnerConstruction::PixelTrk_Barrel(long double r_inn_trkPix, long double r_out_trkPix, G4VisAttributes* VisAtt) 
{
	
	long double l_Pix =  r_out_trkPix/tan(theta_min);
	long double l_Pix_iron =  (r_out_trkPix+widthiron_add)/tan(theta_min);
	

  	G4Tubs *Inner_trkPix = new G4Tubs("Inner_trkPix",r_inn_trkPix,r_out_trkPix, l_Pix, 0,config_obj.max_phi);
	G4LogicalVolume *Inner_trkPix_LV = new G4LogicalVolume(Inner_trkPix, elSi, "Inner_trkPix_LV");
	new G4PVPlacement(
					0,                // no rotation
					G4ThreeVector(0., 0., 0. ),
					Inner_trkPix_LV,          // its logical volume
					"Inner_trkPix_PL",          // its name
					MagFieldLV,                // its mother  volume
					false,            // no boolean operation
					0 ,               // copy number
					fCheckOverlaps
	);
	G4Tubs *Inner_trkPix_iron = new G4Tubs("Inner_trkPix_iron", r_out_trkPix, r_out_trkPix+widthiron_add, l_Pix_iron,0,config_obj.max_phi);
	G4LogicalVolume *Inner_trkPix_LV_iron = new G4LogicalVolume(Inner_trkPix_iron, iron, "Inner_trkPix_LV_iron");
	new G4PVPlacement(
					0,                // no rotation
					G4ThreeVector(0., 0., 0. ),
					Inner_trkPix_LV_iron,          // its logical volume
					"Inner_trkPix_PL_iron",          // its name
					MagFieldLV,                // its mother  volume
					false,            // no boolean operation
					0 ,              // copy number
					fCheckOverlaps
	);
	Inner_trkPix_LV->SetVisAttributes(VisAtt);
	Inner_trkPix_LV_iron->SetVisAttributes(VisAtt);
}

void InnerConstruction::PixelTrk_Barrel(long double r_inn_trkPix, long double r_out_trkPix, long double length_var, G4VisAttributes* VisAtt, bool isOutermostLayer) 
{

	std::string detName = "Inner_trkPix_PL";
	std::string supName = "Inner_trkPix_PL_iron";
	if ( isOutermostLayer ) {
	    detName += "_outermostInner";
	    supName += "_outermostInner";
	}
	
  	G4Tubs *Inner_trkPix = new G4Tubs("Inner_trkPix",r_inn_trkPix,r_out_trkPix, length_var, 0,config_obj.max_phi);
	G4LogicalVolume *Inner_trkPix_LV = new G4LogicalVolume(Inner_trkPix, elSi, "Inner_trkPix_LV");
	new G4PVPlacement(
					0,                // no rotation
					G4ThreeVector(0., 0., 0. ),
					Inner_trkPix_LV,          // its logical volume
					detName,          // its name
					MagFieldLV,                // its mother  volume
					false,            // no boolean operation
					0 ,               // copy number
					fCheckOverlaps
	);
	G4Tubs *Inner_trkPix_iron = new G4Tubs("Inner_trkPix_iron", r_out_trkPix, r_out_trkPix+widthiron_add, length_var,0,config_obj.max_phi);
	G4LogicalVolume *Inner_trkPix_LV_iron = new G4LogicalVolume(Inner_trkPix_iron, iron, "Inner_trkPix_LV_iron");
	new G4PVPlacement(
					0,                // no rotation
					G4ThreeVector(0., 0., 0. ),
					Inner_trkPix_LV_iron,          // its logical volume
					supName,          // its name
					MagFieldLV,               // its mother  volume
					false,            // no boolean operation
					0 ,              // copy number
					fCheckOverlaps
	);
	Inner_trkPix_LV->SetVisAttributes(VisAtt);
	VisAtt->SetForceSolid(true);
	Inner_trkPix_LV_iron->SetVisAttributes(VisAtt);
	
}
