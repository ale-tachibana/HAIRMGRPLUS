# ##### BEGIN GPL LICENSE BLOCK #####
#
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 2
# of the License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software Foundation,
# Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
#
# ##### END GPL LICENCE BLOCK #####

bl_info = {
    "name":"hairmgrplus",
    "author": "Alessandro Tachibana",
    "version": (0,1,0),
    "blender": (2,80,0),
    "location": "Properties",
    "category": "Particle",
    "description": "hair tool",
    "wiki_url": " ",
    "tracker_url":" " 
}

import bpy
import mathutils
import numpy as np
import random


#----------------------------------------
#---------------PROPERTIES---------------
#----------------------------------------
bpy.types.Scene.hamgpImportColl = bpy.props.PointerProperty(
    type = bpy.types.Collection
)
bpy.types.Scene.hamgpExportColl = bpy.props.PointerProperty(
    type = bpy.types.Collection
)




#COPY PASTE FUNCTION VARIABLES
#--------------------------------------------
#hair dynamics
prop_hair_dynamics = [['cloth','settings'],['quality','pin_stiffness']]
prop_hair_dynamics_2 = [[],['use_hair_dynamics']]
prop_hd_structure = [['cloth','settings'],['mass','bending_stiffness','bending_damping']]
prop_hd_structure_2 = [['settings'],['bending_random']]
prop_hd_volume = [['cloth','settings'],['air_damping','air_damping','voxel_cell_size','density_target', 'density_strength']]

#--------------------------------------------
#render
prop_render = [[],['render_type','material_slot','parent']]
prop_render_2 = [['id_data'],['show_instancer_for_render']]
prop_render_path = [['settings'],['use_hair_bspline','render_step']]
prop_render_timing = [['settings'],['use_absolute_path_time','path_start','path_end','length_random']]
prop_render_extra = [['settings'],['use_parent_particles','show_unborn','use_dead']]

#--------------------------------------------
#viewport display
prop_viewport_display = [['settings'] \
                        ,['display_method','display_color','display_percentage' \
                        ,'show_instancer_for_viewport','color_maximum','display_step' \
                        ,'display_size']]

#--------------------------------------------
#children
prop_children = [['settings'] \
                , ['child_type','child_nbr','rendered_child_count' \
                ,'child_length','child_length_threshold','child_seed' \
                ,'child_size', 'child_size_random','child_radius' \
                ,'child_roundness']]

prop_children_parting = [['settings'],['child_parting_factor','child_parting_min','child_parting_max']]
                
prop_children_clumping = [['settings'] \
                        ,['use_clump_curve','clump_factor','clump_shape' \
                        ,'twist','use_twist_curve','twist_curve']]

prop_children_roughness = [['settings'] \
                          ,['use_roughness_curve','roughness_1_size','roughness_endpoint' \
                          ,'roughness_end_shape','roughness_2','roughness_2_size' \
                          , 'roughness_2_threshold', 'roughness_curve']]

prop_children_klin = [['settings'] \
                     ,['kink','kink_amplitude','kink_amplitude_random' \
                     ,'kink_axis','kink_axis_random','kink_frequency' \
                     ,'kink_shape','kink_extra_steps','kink_amplitude' \
                     ,'kink_amplitude_clump','kink_flat','kink_frequency' \
                     ,'kink_shape']]

#--------------------------------------------
prop_hair_shape = [['settings'] \
                  ,['shape','root_radius','tip_radius' \
                  ,'radius_scale','use_close_tip']]


#----------------------------------------
#---------------FUNCTIONS----------------
#----------------------------------------
def debugPrint(string):
    print(string)
    pass

def evaluateObj(name, context):
    despgraph = context.evaluated_depsgraph_get()
    eval_ob = despgraph.objects.get(name, None)
    return eval_ob

def getObjLocal(context):
    return context.scene

def getImpCol(context):
    return getObjLocal(context).hamgpImportColl

def getExpCol(context):
    return getObjLocal(context).hamgpExportColl    

def addToCol(col, obj):
    col.objects.link(obj)
    
def createHairSystem(context, col_name, hair_count):    
    
    hairSys = context.selected_objects[0].modifiers.new(col_name, type='PARTICLE_SYSTEM')
    hairSys.particle_system.settings.type = 'HAIR'
    if hair_count < 4:
        hair_count = 4
    hairSys.particle_system.settings.count = hair_count
    #hairSys.particle_system.is_edited = True
    #hairSys = context.selected_objects[0].particle_systems[-1]
    
    return hairSys

def vect2Array4x1(vector):    
    array = np.array(vector)
    arr4x1 = np.append(array, [1])
    return arr4x1

def calcAngle(da, db):
    da = round(da, 4)
    db = round(db, 4)        
    angle = 0
    if abs(da) > 0:
        angle = np.arctan(abs(db)/abs(da))
        if db > 0 and da > 0: #Q1
            debugPrint('Q1')            
            angle += 0  
        elif db > 0 and da < 0: #Q2
            debugPrint('Q2')
            angle = (np.pi - angle)
        elif db < 0 and da < 0: #Q3
            debugPrint('Q3')           
            angle = -1 * (np.pi - angle)
        elif db < 0 and da > 0: #Q4
            debugPrint('Q4')            
            angle = -1 * (angle)
        elif abs(db) == 0 and da < 0:
            debugPrint('D1')         
            angle = np.pi
    elif db > 0:
        debugPrint('D2')      
        angle = np.pi/2
    elif db < 0:    
        debugPrint('D3')      
        angle = -1 * (np.pi/2) 
    return angle   

def calcRotationMax(position, normal, center):

    debugPrint('----------------------------')
    debugPrint('----------------------------')    
    (dx, dy, dz) = position * -1    
    (nx, ny, nz) = normal
    debugPrint('normal: ' + str(np.round_(normal,4)))    
    ax = 0 #np.arcsin(nx)  #np.arcsin(nx) #* -1d
    ay = 0 #np.arcsin(ny)  #np.deg2rad(-45) #np.arcsin(ny) #* -1
    az = 0 #np.arcsin(nz) * -1

    hh = np.sqrt(np.power(nx,2) + np.power(ny,2))
    
    ax = calcAngle(nz, ny)
    RX = np.array([[1,0,0,0],[0, np.cos(ax), -np.sin(ax), 0],[0, np.sin(ax), np.cos(ax), 0],[0,0,0,1]])
    
    #---------------------------------
    normal_rx = matMul4x4(RX, normal)
    (nrxx, nrxy, nrxz) = normal_rx
    debugPrint('normal_rx: ' + str(np.round_(normal_rx,4)))
    ay = calcAngle(nrxz, nrxx) * -1
    RY = np.array([[np.cos(ay), 0, np.sin(ay),0],[0, 1, 0, 0],[-np.sin(ay), 0, np.cos(ay), 0],[0,0,0,1]])    
    
    #---------------------------------    
    normal_ry = matMul4x4(RY, normal_rx)
    debugPrint('normal_ry: ' + str(np.round_(normal_ry,4)))
    
    
    az = 0 #np.pi
    
    #print('hh: ' + str(hh))             
    debugPrint('----------------------------')
    debugPrint('----------------------------')
    debugPrint('angle x: ' + str(np.rad2deg(ax)))
    debugPrint('angle y: ' + str(np.rad2deg(ay)))
    debugPrint('angle z: ' + str(np.rad2deg(az)))
    
    #ROTATION X - [[1,0,0,0][0, cos Ax, -sin Ax, 0][0, sin Ax, cos Ax, 0][0,0,0,1]]
    #ROTATION Y - [[cos Ay, 0, sin Ay,0][0, 1, 0, 0][-sin Ay, 0, cos Ay, 0][0,0,0,1]]    
    #ROTATION Z - [[cos Az, -sin Az, 0, 0][sin Az, cos Az, 0, 0][0, 0, 1, 0][0,0,0,1]]        
    TM = np.array([[1,0,0,dx],[0,1,0,dy],[0,0,1,dz],[0,0,0,1]])


    RZ = np.array([[np.cos(az), -np.sin(az), 0, 0],[np.sin(az), np.cos(az), 0, 0],[0, 0, 1, 0],[0,0,0,1]])
    
    return TM, RX, RY, RZ
        
def array4x12Vec(array):
    return mathutils.Vector(array[:3])
    
def calcPos(vector, TM, RX, RY, RZ):
    arr4x1 = vect2Array4x1(vector)

    arr4x1 = np.matmul(TM, arr4x1)    
    arr4x1 = np.matmul(RX, arr4x1)    
    arr4x1 = np.matmul(RY, arr4x1)    
    arr4x1 = np.matmul(RZ, arr4x1)
                        
    return array4x12Vec(arr4x1)
        
def printHairKey(context):
    eval_ob = evaluateObj(context.selected_objects[0].name, context)
    evalHair = eval_ob.particle_systems[-1]
    
    debugPrint('----------------------------')
    debugPrint('----------------------------')    
    for index in range(len(evalHair.particles)):        
        hairEval = evalHair.particles[index]
                
        faceindex = eval_ob.closest_point_on_mesh(hairEval.location)[-1]
        position = eval_ob.closest_point_on_mesh(hairEval.location)[1]
        normal = eval_ob.data.polygons[faceindex].normal
        center = eval_ob.data.polygons[faceindex].center
        
        TM, RX, RY, RZ = calcRotationMax(position, normal, center)        
        
        debugPrint('----------------------------')    
        debugPrint('----------------------------') 
        debugPrint('position' + str(position))        
        debugPrint('normal' + str(normal))
        debugPrint('center' + str(eval_ob.data.polygons[faceindex].center))
        
        debugPrint('----------------------------')    
        debugPrint('----------------------------')        
        debugPrint('location: ' + str(hairEval.location))
        debugPrint('velocity: ' + str(hairEval.velocity))      
        debugPrint('rotation: ' + str(hairEval.rotation))        
                  
        for index2 in range(len(hairEval.hair_keys)):
            hairkey = hairEval.hair_keys[index2]
            
            debugPrint('----------------------------')
            debugPrint('co: ' + str(hairkey.co))
            debugPrint('co_local: ' + str(hairkey.co_local))
            
            calc = calcPos(hairkey.co, TM, RX, RY, RZ)
            debugPrint('calc_local: ' + str(calc))    
    
def returnArrayEq(hair_key):
    (xl, yl, zl) = hair_key.co_local
    return xl, yl, zl, hair_key.co

def matMul4x4(AT, co):
    (x, y, z) = co
    pos_array = np.array([x,y,z,1])
            
    co_local = np.matmul(AT,pos_array)
    
    return co_local[:3]


def SetRandomHairPos(hairSys, context):
    
    bpy.ops.particle.particle_edit_toggle()
    bpy.ops.particle.particle_edit_toggle()   
      
    #remove hair with 0 position on any axis
    for index in range(len(hairSys.particles)):
        hair = hairSys.particles[index]
        for index2 in range(1, len(hair.hair_keys)):
            hair_key = hair.hair_keys[index2]
            (x, y, z) = hair_key.co_local
            x = random.uniform(0.01,0.03) * index2           
            y = random.uniform(0.1,0.3)  * index2             
            z = random.uniform(-0.1,-0.3) * index2
            hair_key.co_local =  mathutils.Vector([x, y, z])
            

def Curves2Hair(collCurves, hairSys, context):        
    #hairSys = context.selected_objects[0].particle_systems[-1]        
    bpy.ops.particle.particle_edit_toggle()
    bpy.ops.particle.particle_edit_toggle() 
            
    print(hairSys.name)
    eval_ob = evaluateObj(context.selected_objects[0].name, context)
    evalHair = eval_ob.particle_systems[-1]
    
    for index in range(len(hairSys.particles)):
        #print(index)
        if index <= len(collCurves):
            curve = collCurves[index]
                        
            #curveEval = evaluateObj(curve.name, context)
            hair = hairSys.particles[index]
            hairEval = evalHair.particles[index]
                        
            bezier_pts = curve.data.splines.items()[0][1].bezier_points.items()            
            
            
            faceindex = eval_ob.closest_point_on_mesh(hairEval.location)[-1]
            position = eval_ob.closest_point_on_mesh(hairEval.location)[1]
            normal = eval_ob.data.polygons[faceindex].normal
            center = eval_ob.data.polygons[faceindex].center
                                    
            #-------------------------------------
            #calculate the transformation matrix            
            #-------------------------------------                        
            TM, RX, RY, RZ = calcRotationMax(position, normal, center) 
                                    
            for index2 in range(len(hair.hair_keys)):
                if index2 <= len(bezier_pts) - 1:
                    points = bezier_pts[index2][1]
                    
                    hair_key = hair.hair_keys[index2]
                    hair_kEval = hairEval.hair_keys[index2]
                                        
                    calc = calcPos(points.co, TM, RX, RY, RZ)
                    
                    #FUCK THIS MOTHERFUCKING SHIT
                    #WHY CAN'T YOU JUST SET THE .CO VARIABLE????????
                    #FUCK FUCK FUCK!!!!!!!!!!!!!!!!!!
                    #hair_key.co_local = calc
                    hair_key.co_local = calc


    
def Hair2Curves(hairSys, collCurves, context):
    for hair in hairSys.particles.items():        
        #create curve
        newCurve = bpy.data.curves.new('haircurve', type='CURVE')
        newCurve.dimensions = '3D'
                
        spline = newCurve.splines.new('BEZIER')
        spline.bezier_points.add(len(hair[1].hair_keys.items()) -1)
        #print('---------------------------------------------') 
        
        for index, hairkey in hair[1].hair_keys.items():
            #change curves
            spline.bezier_points[index].co = hairkey.co
            spline.bezier_points[index].handle_left_type = "AUTO"
            spline.bezier_points[index].handle_right_type = "AUTO"
            #print(hairkey[1].co_local)
            #print((x, y, z))
                                        
        curveObj = bpy.data.objects.new('curveobj', newCurve)
        #curveObj.rotation_mode = 'XYZ'
        #curveObj.rotation_quaternion = hair[1].rotation
        #curveObj.location = hair[1].location                
                
        addToCol(collCurves, curveObj)

#----------------------------------------
#----------------------------------------
#-------------OBJ TYPE-------------------
#----------------------------------------
#----------------------------------------



#----------------------------------------
#----------------------------------------
#---------------INTERFACE----------------
#----------------------------------------
#----------------------------------------
class hairmgrplusPanel:    
    bl_space_type = "PROPERTIES"
    bl_region_type = "WINDOW"
    bl_context = "particle" 
    bl_options = {'DEFAULT_CLOSED'}
    

#----------------------------------------
#------------OT--------------------------
#----------------------------------------
class HAIRMGRPLUS_OT_force_enable_adv(bpy.types.Operator):
    bl_idname = "hairmgrplus.force_enable_adv"    
    bl_label = "TOGGLE ENABLE ADVANCED HAIR"
    
    def execute(self, context):
        if len(context.selected_objects) == 0:
            return {'CANCELLED'}
        elif context.selected_objects[0].particle_systems is None:
            return {'CANCELLED'}
        elif context.selected_objects[0].particle_systems.active is None:
            return {'CANCELLED'}
        elif context.selected_objects[0].particle_systems.active.settings.type != 'HAIR':
            return {'CANCELLED'}
        
        settings = context.selected_objects[0].particle_systems.active.settings
        settings.use_advanced_hair = not settings.use_advanced_hair
        
        return {'FINISHED'}

        
class HAIRMGRPLUS_OT_copy_parameters(bpy.types.Operator):
    bl_idname = "hairmgrplus.copy_parameters"    
    bl_label = "COPY PARAMETERS"   
    
    @classmethod
    def poll(cls, context):
        return False
    
    def execute(self, context):
        return {'FINISHED'}    
    
    
class HAIRMGRPLUS_OT_paste_parameters(bpy.types.Operator):
    bl_idname = "hairmgrplus.paste_parameters"    
    bl_label = "PASTE PARAMETERS"   
    
    @classmethod
    def poll(cls, context):
        return False
    
    def execute(self, context):
        return {'FINISHED'}        
    
class HAIRMGRPLUS_OT_create_hair_sys(bpy.types.Operator):
    bl_idname = "hairmgrplus.create_hair_sys"    
    bl_label = "Create Hair System"
    
    @classmethod
    def poll(cls, context):
        var_return = True
        if len(context.selected_objects) == 0:
            var_return = False    
        elif getImpCol(context) is None:
            var_return = False
        elif len(getImpCol(context).all_objects.items()) == 0:
            var_return = False
            
        return var_return
    
    def execute(self, context): 
        self.obj = context.object
                
        importCol = getImpCol(context)
        filteredCol = []
        
        for obj in importCol.all_objects.items():
            if obj[1].type == 'CURVE':
                filteredCol.append(obj[1])
        
        if len(filteredCol) == 0:
            print('Objects need to be of the type = CURVE')
            return {'CANCELLED'}
        
        pointsCount = 0
        for obj in filteredCol:
            bezier_pts = obj.data.splines.items()[0][1].bezier_points.items()
            if len(bezier_pts) > pointsCount:
                pointsCount = len(bezier_pts)                
        
        partSystem = createHairSystem(context, importCol.name, len(filteredCol))
        
        partSystem.particle_system.settings.hair_step = pointsCount - 1
                
        return {'FINISHED'}

class HAIRMGRPLUS_OT_test_calc(bpy.types.Operator):        
    bl_idname = "hairmgrplus.test_calc"    
    bl_label = "Test Calc"
    
    @classmethod
    def poll(cls, context):
        var_return = True
        if len(context.selected_objects) == 0:
            var_return = False    
        elif getImpCol(context) is None:
            var_return = False
        elif len(getImpCol(context).all_objects.items()) == 0:
            var_return = False
            
        return var_return
        
    def execute(self, context): 
        self.obj = context.object
        
        hairSys = context.selected_objects[0].particle_systems[-1]
        
        #print(dir(hairSys))
        #SetRandomHairPos(hairSys, context)  
        printHairKey(context)
        #ExportData(hairSys, context)              
        return {'FINISHED'}
    
    
class HAIRMGRPLUS_OT_import_from_curves(bpy.types.Operator):
    bl_idname = "hairmgrplus.import_from_curves"    
    bl_label = "Import from Curves"
    
    @classmethod
    def poll(cls, context):
        var_return = True
        if len(context.selected_objects) == 0:
            var_return = False    
        elif getImpCol(context) is None:
            var_return = False
        elif len(getImpCol(context).all_objects.items()) == 0:
            var_return = False
            
        return var_return
    
    def execute(self, context): 
        self.obj = context.object
                
        importCol = getImpCol(context)
        filteredCol = []
        
        for obj in importCol.all_objects.items():
            if obj[1].type == 'CURVE':
                filteredCol.append(obj[1])
        
        if len(filteredCol) == 0:
            print('Objects need to be of the type = CURVE')
            return {'CANCELLED'}        
        
        hairSys = context.selected_objects[0].particle_systems[-1]
        
        #print(dir(hairSys))    
        Curves2Hair(filteredCol, hairSys, context)
                
        return {'FINISHED'}


class HAIRMGRPLUS_OT_export_to_curves(bpy.types.Operator):    
    bl_idname = "hairmgrplus.export_to_curves"
    bl_label = "Export to Curves"

    @classmethod
    def poll(cls, context):
        var_return = True
        if len(context.selected_objects) == 0:
            var_return = False            
        elif context.selected_objects[0].particle_systems is None:
            var_return = False
        elif context.selected_objects[0].particle_systems.active is None:
            var_return = False 
        elif context.selected_objects[0].particle_systems.active.settings.type != 'HAIR':
            var_return = False
        elif getExpCol(context) is None:
            var_return = False 
        
        return var_return
        
    def execute(self, context):
        self.obj = context.object
        
        bpy.ops.particle.particle_edit_toggle()
        bpy.ops.particle.particle_edit_toggle()
                
        if len(context.selected_objects[0].particle_systems.active.particles.items()) == 0:
            print('Particle system has no hairs')
            return {'CANCELLED'}
        
        #the .co fields on the hair_keys won't be loaded unless evaluated before
        #--------------------------------------------
        eval_ob = evaluateObj(context.selected_objects[0].name, context)
        #--------------------------------------------
                
        Hair2Curves(eval_ob.particle_systems.active, getExpCol(context), context)
        
        return {'FINISHED'}



#----------------------------------------
#------------PT--------------------------
#----------------------------------------

class HAIRMGRPLUS_PT_panel(hairmgrplusPanel, bpy.types.Panel): 
    bl_label = "Hair Manager Plus"   
    #bl_owner_id = "HAIRMGRPLUS_PT_panel"
    bl_options = {'DEFAULT_CLOSED'}
    
    def draw(self, context):
        layout = self.layout
            
        row = layout.row()
        #row.label(text = "-------")
        
        
class HAIRMGRPLUS_PT_manager(hairmgrplusPanel, bpy.types.Panel):
    bl_label = "Hair Manager"
    bl_parent_id = "HAIRMGRPLUS_PT_panel"
        
    def draw(self, context):
        layout = self.layout
        
        box = layout.box()
        
        row = box.row()        
        row.label(text = "STILL UNDER WORK")
                
        row = layout.row()
        row.operator("hairmgrplus.force_enable_adv") 
        
        row = layout.row()
        row.operator("hairmgrplus.copy_parameters", text="-----------") 
        
        row = layout.row()
        row.operator("hairmgrplus.paste_parameters", text="------------") 


class HAIRMGRPLUS_PT_import_curves(hairmgrplusPanel, bpy.types.Panel):
    bl_label = "Import Curves"
    bl_parent_id = "HAIRMGRPLUS_PT_panel"
            
    def draw(self, context):
        self.obj = context.object
        layout = self.layout
                
        row = layout.row()        
        row.label(text = "DOES NOT WORK")        
        
        row = layout.row()
        row.prop(getObjLocal(context), 'hamgpImportColl', text = "Import Col.")

        
        row = layout.row()
        row.operator("hairmgrplus.create_hair_sys", text="Create Hair")         
                
        row = layout.row()
        row.operator("hairmgrplus.test_calc", text="Test Transf. Matrix")                
                
        row = layout.row()
        row.operator("hairmgrplus.import_from_curves", text="Import from Curves")         

    
class HAIRMGRPLUS_PT_export_curves(hairmgrplusPanel, bpy.types.Panel):
    bl_label = "Export Curves"
    bl_parent_id = 'HAIRMGRPLUS_PT_panel'
        
    def draw(self, context):
        self.obj = context.object
        layout = self.layout
                
        row = layout.row()        
        row.label(text = "KINDA OFF WORK")
        
        #box = layout.box()            
        row = layout.row()        
        row.prop(getObjLocal(context), 'hamgpExportColl', text = "Export Col.")

        row = layout.row()        
        row.operator("hairmgrplus.export_to_curves", text="Export to Curves")       

#----------------------------------------
#----------------------------------------
#----------------------------------------




#----------------------------------------
#----------------REGISTER----------------
#----------------------------------------
classes = (    
    HAIRMGRPLUS_PT_panel
    , HAIRMGRPLUS_OT_force_enable_adv
    , HAIRMGRPLUS_OT_copy_parameters
    , HAIRMGRPLUS_OT_paste_parameters
    , HAIRMGRPLUS_OT_create_hair_sys
    , HAIRMGRPLUS_OT_test_calc
    , HAIRMGRPLUS_OT_import_from_curves
    , HAIRMGRPLUS_OT_export_to_curves    
    , HAIRMGRPLUS_PT_manager
    , HAIRMGRPLUS_PT_import_curves
    , HAIRMGRPLUS_PT_export_curves    
)

def register():        
    for cls in classes:
        bpy.utils.register_class(cls) 
        
def unregister():        
    for cls in reversed(classes):
        bpy.utils.unregister_class(cls)

if __name__ == '__main__':
    register()
    