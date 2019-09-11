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
import ast

#----------------------------------------
#---------------PROPERTIES---------------
#----------------------------------------
bpy.types.Scene.hamgpImportColl = bpy.props.PointerProperty(
    type = bpy.types.Collection
)
bpy.types.Scene.hamgpExportColl = bpy.props.PointerProperty(
    type = bpy.types.Collection
)


#----------------------------------------
#----------------------------------------
#COPY PASTE FUNCTION VARIABLES
#----------------------------------------
#----------------------------------------
COPY_PASTE_KEY = 'HMGRP_COPY_DATA'

#hair dynamics
prop_hair_dynamics = [[],['use_hair_dynamics']]
prop_hair_dynamics_2 = [['cloth','settings'],['quality','pin_stiffness']]

prop_hd_structure = [['cloth','settings'],['mass','bending_stiffness','bending_damping']]
prop_hd_structure_2 = [['settings'],['bending_random']]
prop_hd_volume = [['cloth','settings'],['air_damping','air_damping','voxel_cell_size','density_target', 'density_strength','internal_friction']]

#--------------------------------------------
#render
prop_render = [['settings'],['render_type','material_slot']]
prop_render_2 = [['id_data'],['show_instancer_for_render']]
prop_render_3 = [[],['parent']]
prop_render_path = [['settings'],['use_hair_bspline','render_step']]
prop_render_timing = [['settings'],['use_absolute_path_time','path_start','path_end','length_random']]
prop_render_extra = [['settings'],['use_parent_particles','show_unborn','use_dead']]

#--------------------------------------------
#viewport display
prop_viewport_display = [['settings'] \
                        ,['display_method','display_color','display_percentage' \
                        ,'display_size','color_maximum','display_step' ]]
                        
prop_viewport_display_2 = [['id_data'],['show_instancer_for_viewport']]                      
                        
#--------------------------------------------
#children
prop_children = [['settings'] \
                , ['child_type','child_nbr','rendered_child_count' \
                ,'child_length','child_length_threshold' \
                ,'child_size', 'child_size_random','child_radius' \
                ,'child_roundness']]
                
prop_children_2 = [[],['child_seed']]                

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
#shape
prop_hair_shape = [['settings'] \
                  ,['shape','root_radius','tip_radius' \
                  ,'radius_scale','use_close_tip']]
                  
#---------------------------------------------
#velocity
prop_velocity = [['settings'] \
                ,['normal_factor','tangent_factor', 'tangent_phase' \
                ,'object_align_factor','particle_factor','object_factor' \
                , 'factor_random']]                  

#----------------------------------------
#----------------------------------------
#add menu name
prop_all_hd = ['Hair Dynamics',[prop_hair_dynamics, prop_hair_dynamics_2, prop_hd_structure, prop_hd_structure_2, prop_hd_volume]]
prop_all_render = ['Render',[prop_render, prop_render_2, prop_render_3, prop_render_path, prop_render_timing, prop_render_extra]]
prop_all_vd = ['Viewport Display',[prop_viewport_display, prop_viewport_display_2]]
prop_all_children = ['Children',[prop_children, prop_children_parting, prop_children_clumping, prop_children_roughness, prop_children_klin]]
prop_all_hs = ['Hair Shape',[prop_hair_shape]]
prop_all_veloc = ['Velocity',[prop_velocity]]
#prop_test_child = ['Test', [prop_children_roughness]]

#prop_all is used to build the panel COPY buttons
prop_all = [prop_all_hd, prop_all_render, prop_all_vd, prop_all_children, prop_all_hs, prop_all_veloc]



#----------------------------------------
#---------------FUNCTIONS----------------
#----------------------------------------
def debugPrint(string):
    print(string)
    pass

def errorPrint(string):
    print(string)
    pass

def setClipBoard(var):
    bpy.context.window_manager.clipboard = var

def getClipBoard():
    return bpy.context.window_manager.clipboard

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
    
def IsParticleSelected(context):
    var_return = True
    if len(context.selected_objects) == 0:
        var_return = False            
    elif context.selected_objects[0].particle_systems is None:
        var_return = False
    elif context.selected_objects[0].particle_systems.active is None:
        var_return = False 
    elif context.selected_objects[0].particle_systems.active.settings.type != 'HAIR':
        var_return = False
    return var_return

#----------------------------------------
#----------------------------------------    
#----------------------------------------
class Import:    
    @classmethod    
    def createHairSystem(self, context, col_name, hair_count):    
        
        hairSys = context.selected_objects[0].modifiers.new(col_name, type='PARTICLE_SYSTEM')
        hairSys.particle_system.settings.type = 'HAIR'
        if hair_count < 4:
            hair_count = 4
        hairSys.particle_system.settings.count = hair_count
        #hairSys.particle_system.is_edited = True
        #hairSys = context.selected_objects[0].particle_systems[-1]
        
        return hairSys
        
    @classmethod    
    def __calcAngle(self, da, db):
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
    
    @classmethod
    def __solveRotationMax(self, position, particle):        
        A = B = C = D = E = F = G = H = I = J = K = 0
        
        (dx, dy, dz) = position * -1
        TM = np.array([[1,0,0,dx],[0,1,0,dy],[0,0,1,dz],[0,0,0,1]])
        
        co1 = self.__matMul4x4(TM, particle.hair_keys[1].co)
        co2 = self.__matMul4x4(TM, particle.hair_keys[2].co)
        co3 = self.__matMul4x4(TM, particle.hair_keys[3].co)                
                
        (x1, y1, z1) = co1
        (x2, y2, z2) = co2
        (x3, y3, z3) = co3
        (xl1, yl1, zl1) = particle.hair_keys[1].co_local
        (xl2, yl2, zl2) = particle.hair_keys[2].co_local
        (xl3, yl3, zl3) = particle.hair_keys[3].co_local
        
        a1 = np.array([[x1, y1, z1],[x2, y2, z2], [x3, y3, z3]])
        b1 = np.array([xl1, xl2, xl3])
        
        a2 = np.array([[x1, y1, z1],[x2, y2, z2], [x3, y3, z3]])
        b2 = np.array([yl1, yl2, yl3])
        
        a3 = np.array([[x1, y1, z1],[x2, y2, z2], [x3, y3, z3]])
        b3 = np.array([zl1, zl2, zl3])

        try:
            r1 = np.linalg.solve(a1,b1)
            (A, B, C) = r1
        except Exception as e:
            errorPrint(str(e))        
        
        try:
            r2 = np.linalg.solve(a2,b2)
            (D, E, F) = r2
        except Exception as e:
            errorPrint(str(e))                
        
        try:
            r3 = np.linalg.solve(a3,b3)
            (H, I, J) = r3
        except Exception as e:
            errorPrint(str(e))                
                
        MATRIX = [[A, B, C, dx],[D, E, F, dy],[H, I, J, dz],[0, 0, 0, 1]]
        
        return MATRIX
    
    @classmethod
    def __findZMatRot(self, particle, TM, RX, RY):
        #|cosA  sinA   0    0 
        #|sinA  cosA   0    0
        #|0      0     1    0
        #|0      0     0    1
        RZ = [[1, 0, 0, 0],[0, 1, 0, 0],[0, 0, 1, 0],[1, 0, 0, 1]]
        pos1 = self.__calcPos(particle.hair_keys[1].co, TM, RX, RY, RZ)
        pos2 = self.__calcPos(particle.hair_keys[2].co, TM, RX, RY, RZ)
        
        az = 0
        cosZ = 0
        sinZ = 0
        (x1, y1, z1) = pos1
        (x2, y2, z2) = pos2
        
        (xl1, yl1, zl1) = particle.hair_keys[1].co_local
        (xl2, yl2, zl2) = particle.hair_keys[2].co_local
        
        # x_local = cozZ * x -sinZ * y
        a = np.array([[x1, y1],[x2, y2]])
        b = np.array([xl1, xl2])
        try:
            x = np.linalg.solve(a,b)
            (cosZ, sinZ) = x                  
        except Exception as e:
            errorPrint(str(e))
        
        try:
            az1 = np.arccos(cosZ) 
            if not np.isnan(az1):
                az = az1
        except Exception as e:
            errorPrint(str(e))
        
        try:
            az2 = np.arcsin(sinZ)
            if not np.isnan(az2):
                az = az2 
        except Exception as e:
            errorPrint(str(e))
        
        RZ = np.array([[np.cos(az), -np.sin(az), 0, 0],[np.sin(az), np.cos(az), 0, 0],[0, 0, 1, 0],[0,0,0,1]])
        debugPrint('z angle: ' + str(az))   
        debugPrint('RZ: ' + str(RZ))
        return RZ
    
    @classmethod    
    def __calcRotationMax(self, position, normal, center):

        debugPrint('----------------------------')
        debugPrint('----------------------------')      
        (dx, dy, dz) = position * -1    
        (nx, ny, nz) = normal
        debugPrint('normal: ' + str(np.round_(normal,4)))    
        ax = 0 #np.arcsin(nx)  #np.arcsin(nx) #* -1d
        ay = 0 #np.arcsin(ny)  #np.deg2rad(-45) #np.arcsin(ny) #* -1
        az = 0 #np.arcsin(nz) * -1
        
        hh = np.sqrt(np.power(nx,2) + np.power(ny,2))
        
        ax = self.__calcAngle(nz, ny)
        RX = np.array([[1,0,0,0],[0, np.cos(ax), -np.sin(ax), 0],[0, np.sin(ax), np.cos(ax), 0],[0,0,0,1]])
        
        #---------------------------------
        normal_rx = self.__matMul4x4(RX, normal)
        (nrxx, nrxy, nrxz) = normal_rx
        debugPrint('normal_rx: ' + str(np.round_(normal_rx,4)))
        ay = self.__calcAngle(nrxz, nrxx) * -1
        RY = np.array([[np.cos(ay), 0, np.sin(ay),0],[0, 1, 0, 0],[-np.sin(ay), 0, np.cos(ay), 0],[0,0,0,1]])    
        
        #---------------------------------    
        normal_ry = self.__matMul4x4(RY, normal_rx)
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
    
    @classmethod
    def __vect2Array4x1(self, vector):    
        array = np.array(vector)
        arr4x1 = np.append(array, [1])
        return arr4x1    
 
    @classmethod
    def __array4x12Vec(self, array):
        return mathutils.Vector(array[:3])

    @classmethod
    def __calcPos(self, vector, TM, RX, RY, RZ):
        arr4x1 = self.__vect2Array4x1(vector)

        arr4x1 = np.matmul(TM, arr4x1)    
        arr4x1 = np.matmul(RX, arr4x1)    
        arr4x1 = np.matmul(RY, arr4x1)    
        arr4x1 = np.matmul(RZ, arr4x1)
                            
        return self.__array4x12Vec(arr4x1)

    @classmethod                
    def printHairKey(self, context):
        eval_ob = evaluateObj(context.selected_objects[0].name, context)
        evalHair = eval_ob.particle_systems[-1]
        
        debugPrint('----------------------------')
        debugPrint('----------------------------')    
        for index in range(len(evalHair.particles)):        
            hairEval = evalHair.particles[index]
                    
            faceindex = eval_ob.closest_point_on_mesh(hairEval.location)[-1]
            position = eval_ob.closest_point_on_mesh(hairEval.location)[1]            
            #normal = eval_ob.data.polygons[faceindex].normal
            normal = hairEval.velocity
            center = eval_ob.data.polygons[faceindex].center
            
            #TM, RX, RY, RZ = self.__calcRotationMax(position, normal, center)        
            #RZ = self.__findZMatRot(hairEval, TM, RX, RY)
            
            TMAT = self.__solveRotationMax(position, hairEval)             
            
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
                
                #calc = self.__calcPos(hairkey.co, TM, RX, RY, RZ)
                calc = self.__matMul4x4(TMAT, hairkey.co)
                
                debugPrint('calc_local: ' + str(calc))    

    @classmethod            
    def returnArrayEq(self, hair_key):
        (xl, yl, zl) = hair_key.co_local
        return xl, yl, zl, hair_key.co

    @classmethod    
    def __matMul4x4(self, AT, co):
        (x, y, z) = co
        pos_array = np.array([x,y,z,1])
                
        co_local = np.matmul(AT,pos_array)
        
        return co_local[:3]

    @classmethod    
    def SetRandomHairPos(self, hairSys, context):
        
        bpy.ops.particle.particle_edit_toggle()
        bpy.ops.particle.particle_edit_toggle()   
          
        #remove hair with 0 position on any axis
        for index in range(len(hairSys.particles)):
            hair = hairSys.particles[index]
                        
            for index2 in range(1, len(hair.hair_keys)):
                hair_key = hair.hair_keys[index2]
                (x, y, z) = hair_key.co_local
                x = random.uniform(-0.3,0.3) * index2           
                y = random.uniform(-0.3,0.3) * index2             
                z = random.uniform(0.1,0.3) * index2
                hair_key.co_local =  mathutils.Vector([x, y, z])
            
            hair.hair_keys[0].co_local = mathutils.Vector([0, 0, 0])
                
    @classmethod    
    def Curves2Hair(self, collCurves, hairSys, context):        
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
                #position = eval_ob.closest_point_on_mesh(hairEval.location)[1]
                position = hairEval.hair_keys[0].co #uses the first hair key as the translation delta, works if co_local = 0
                normal = eval_ob.data.polygons[faceindex].normal
                center = eval_ob.data.polygons[faceindex].center
                
                #-------------------------------------
                #calculate the transformation matrix            
                #-------------------------------------                        
                #TM, RX, RY, RZ = self.__calcRotationMax(position, normal, center)               
                #RZ = self.__findZMatRot(hair, TM, RX, RY)
                
                TMAT = self.__solveRotationMax(position, hair)  
                
                for index2 in range(len(hair.hair_keys)):
                    if index2 <= len(bezier_pts) - 1:
                        points = bezier_pts[index2][1]
                        
                        hair_key = hair.hair_keys[index2]
                        hair_kEval = hairEval.hair_keys[index2]
                                            
                        #calc = self.__calcPos(points.co, TM, RX, RY, RZ)
                        calc = self.__matMul4x4(TMAT, points.co)
                        
                        #FUCK THIS MOTHERFUCKING SHIT
                        #WHY CAN'T YOU JUST SET THE .CO VARIABLE????????
                        #FUCK FUCK FUCK!!!!!!!!!!!!!!!!!!
                        #hair_key.co_local = calc
                        hair_key.co_local = calc


#----------------------------------------
#----------------------------------------    
#----------------------------------------
class Export:
    @classmethod            
    def Hair2Curves(self, hairSys, collCurves, context):
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
            
    @classmethod
    def Hair2Poly(self, hairSys, collCurves, context):
        pass

#----------------------------------------
#---------------------------------------



#----------------------------------------
#----------------------------------------
#-------------OBJ TYPE-------------------
#----------------------------------------
#----------------------------------------
class HMGRPCopyPasteType:
        
    class CopyPasteData:
        
        def __init__(self, path, property, value, vartype, children):            
            self.path = path
            self.property = property
            self.value = value
            self.vartype = vartype
            self.children = children
                                            
        def toArray(self):
            vchildren =[]
            for ch in self.children:
                vchildren.append(ch.toArray())
            
            return [self.path, self.property, self.value, self.vartype, vchildren]
        
        @classmethod
        def arrayToString(self, items):
            returnstr = '['
            for i in range(len(items)): 
                if i > 0:
                    returnstr += ','
                returnstr += "\"" + str(items[i]) + "\""
            returnstr += ']'
            return returnstr
        
        @classmethod
        def castString(self, value, typedef):
            if str(typedef) == '<class \'NoneType\'>':
                return None
            elif str(typedef) == '<class \'bool\'>':
                if value == 'True':
                    return True
                else:
                    return False
            elif str(typedef) == '<class \'float\'>':
                return float(value)
            elif str(typedef) == '<class \'int\'>':
                return int(value)
            elif str(typedef) == '<class \'Vector\'>':
                tmpstr = value[value.find("(")+1:value.find(")")]
                tmparr = tmpstr.split(', ')
                (x, y, z) = tmparr
                x = float(x)
                y = float(y)
                z = float(z)                               
                tmpvec = mathutils.Vector([x,y,z])
                #debugPrint(tmpvec)
                return mathutils.Vector(tmpvec)
            else:        
                return value
    
    #----------------------------------------------------
    #----------------------------------------------------
    def __CPDObjFromArray(self, parseArray):
        vpath = ''
        vproperty = ''
        vvalue = ''
        vvartype = ''
        vchildren = []
        
        vtmpch = []
        debugPrint('parseArray: ' + str(parseArray))       
        if len(parseArray) == 5:
            vpath = parseArray[0]
            vproperty = parseArray[1]
            vvalue = parseArray[2]
            vvartype = parseArray[3]
            try:
                vtmpch = ast.literal_eval(parseArray[4])
            except:
                pass
        
        if str(vvartype) == '<class \'bpy.types.CurveMapping\'>' and len(vtmpch) > 0:
            for item in vtmpch:
                debugPrint('------------------------' + item)
                vchild = self.__CPDObjFromArray(item)
                vchildren.append(vchild)
                
            
        return self.CopyPasteData(vpath, vproperty, vvalue, vvartype, vchildren)
    
    def __CPDObjFromValues(self, path, property, value):
        vpath = str(path)
        vproperty = str(property)
        vvalue = str(value)
        vvartype = type(value)
        vchildren = []        
        
        '''
        if str(vvartype) == '<class \'bpy.types.CurveMapping\'>' or str(vvartype) == '<class \'bpy_prop_collection\'>':
            tmp_path = path.copy()
            tmp_path.append(vproperty)
            
            #debugPrint(tmp_path)
            #debugPrint(dir(value))            
            for tmpprop in dir(value):
                if str(tmpprop) != 'bl_rna' and str(tmpprop) != 'rna_type':
                    tmpvalue = getattr(value, tmpprop)
                    #debugPrint(dir(value))                            
                    data = self.__CPDObjFromValues(tmp_path, tmpprop, tmpvalue)
                    #debugPrint('data: ' + str(data.toArray()))   
                    vchildren.append(data)
        '''
        return self.CopyPasteData(vpath, vproperty, vvalue, vvartype, vchildren)
        
    
    def __init__(self, strData=''):
        #debugPrint(strData)
        self.__arrayCPData = []
        if len(strData) > 0:
            #debugPrint('paste data: ' + str(strData))
            try:
                tmpData = ast.literal_eval(strData)
                #print(tmpData[0])
                if tmpData[0] == COPY_PASTE_KEY:
                    #debugPrint(tmpData[1])
                    self.__arrayCPData = tmpData[1]
            except Exception as e:
                errorPrint('error:' + str(e))
    
    def __loadChildren(self, childrenArr):
        return childrenArr
    
    def __procPath(self, var, hairSys):
        returnVar = hairSys
        for path in var:
            try:
                returnVar = getattr(returnVar, path)
            except:
                errorPrint('error path: ' + str(var))
        return returnVar
    
    def Load(self, dataArr, hairSys):        
        for data in dataArr:
            objpath = self.__procPath(data[0], hairSys)            
            for prop in data[1]:                
                value = None                                
                try:
                    value = getattr(objpath, prop)                               
                except Exception as e:
                    errorPrint('error:' + str(e))
                
                cpdata = self.__CPDObjFromValues(data[0], prop, value)                                        
                cpstring = cpdata.toArray() 
                self.__arrayCPData.append(cpstring)                             
    
    def __setHairSysParam(self, hairSys, cpdata):
        #cpdata = CopyPasteData
        #debugPrint('path: ' + cpdata.path) 
        if cpdata.property != '':       
            objpath = self.__procPath(ast.literal_eval(cpdata.path), hairSys)
            #debugPrint(objpath)
            try:
                #setattr(object, name, value)
                value = self.CopyPasteData.castString(cpdata.value, cpdata.vartype)
                setattr(objpath, cpdata.property, value)
                #debugPrint(str(cpdata.property) + ':' + str(value))
            except Exception as e:
                errorPrint('error:' + str(e) + ', parameter: ' + str(cpdata.toArray()))
            
    def __UpdtItem2Hair(self, hairSys, arrayData):
        #debugPrint(arrayData)
        for item in arrayData:
            #debugPrint(item)
            cpdata = self.__CPDObjFromArray(item)
            self.__setHairSysParam(hairSys, cpdata)    
            if len(cpdata.children) > 0:
                self.__UpdtItem2Hair(hairSys, cpdata.children)
    
    def UpdateHair(self, hairSys):
        self.__UpdtItem2Hair(hairSys, self.__arrayCPData)
    
    def Print(self):                
        for item in self.__arrayCPData:
            debugPrint(str(item))
            debugPrint('-------------------------------')
    
    def ToString(self):
        string =  '[\"' + COPY_PASTE_KEY + '\",['
        for i in range(len(self.__arrayCPData)):
            if i > 0:
                string += ','
            string += self.CopyPasteData.arrayToString(self.__arrayCPData[i])
        string += ']]'
        return string
    #----------------------------------------------------
    #----------------------------------------------------
    
        
    
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
    
    @classmethod
    def poll(cls, context):
        return IsParticleSelected(context)
        
    def execute(self, context):                    
        settings = context.selected_objects[0].particle_systems.active.settings
        settings.use_advanced_hair = not settings.use_advanced_hair
        
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
        
        partSystem = Import.createHairSystem(context, importCol.name, len(filteredCol))
        
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
        Import.printHairKey(context)
        #ExportData(hairSys, context)              
        return {'FINISHED'}


class HAIRMGRPLUS_OT_random_pos(bpy.types.Operator):        
    bl_idname = "hairmgrplus.random_pos"    
    bl_label = "Random Pos"
    
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
        Import.SetRandomHairPos(hairSys, context)  
        #Import.printHairKey(context)
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
        Import.Curves2Hair(filteredCol, hairSys, context)
                
        return {'FINISHED'}


class HAIRMGRPLUS_OT_export(bpy.types.Operator):    
    bl_idname = "hairmgrplus.export"
    bl_label = "Export"
    
    export_type: bpy.props.StringProperty()    

    @classmethod
    def poll(cls, context):
        var_return = True
        if not IsParticleSelected(context):
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
        
        if self.export_type == 'CURVE': 
            Export.Hair2Curves(eval_ob.particle_systems.active, getExpCol(context), context)
        elif self.export_type == 'POLY':
            Export.Hair2Poly(eval_ob.particle_systems.active, getExpCol(context), context)
            
        
        return {'FINISHED'}


#----------------------------------------a
#--------COPY/PASTE PARAMETERS-----------
#----------------------------------------
class HAIRMGRPLUS_OT_copy_parameters(bpy.types.Operator):
    bl_idname = "hairmgrplus.copy_parameters"    
    bl_label = "COPY PARAMETERS"

    parameter_copy: bpy.props.StringProperty()
            
    @classmethod
    def poll(cls, context):
        return IsParticleSelected(context)        
        
    def execute(self, context):    
        #debugPrint(self.parameter_copy)
        hairSys = context.selected_objects[0].particle_systems.active
        
        copyData = HMGRPCopyPasteType()
                
        if self.parameter_copy == 'ALL':
            all = []
            
            for pp in prop_all:
                for pp2 in pp[1]:
                    all.append(pp2)
                                                   
            copyData.Load(all, hairSys)
        else:
            index = int(self.parameter_copy)            
            copyData.Load(prop_all[index][1], hairSys)            
        #copyData.Print()
        setClipBoard(copyData.ToString())
        
        return {'FINISHED'}

class HAIRMGRPLUS_OT_paste_parameters(bpy.types.Operator):
    bl_idname = "hairmgrplus.paste_parameters"
    bl_label = "PASTE PARAMETERS"    
    
    @classmethod
    def poll(cls, context):
        return IsParticleSelected(context)
        
    def execute(self, context):
        data = getClipBoard()
        hairSys = context.selected_objects[0].particle_systems.active
        
        copyData = HMGRPCopyPasteType(data)
        copyData.UpdateHair(hairSys)
        
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
        
        
class HAIRMGRPLUS_PT_manager(hairmgrplusPanel, bpy.types.Panel):
    bl_label =  "Management Panel"
    bl_parent_id = "HAIRMGRPLUS_PT_panel"
    bl_options = {'DEFAULT_CLOSED'}
    
    def draw(self, context):
        layout = self.layout
        
        #box = layout.box()
            
        row = layout.row()
        row.operator("hairmgrplus.force_enable_adv")         
    
        layout.separator()    
        layout.separator()        
                
        row = layout.row()
        row.operator("hairmgrplus.paste_parameters", text="PASTE PARAMETERS") 

        layout.separator()
        
        row = layout.row()
        row.operator("hairmgrplus.copy_parameters", text="COPY PARAMETERS").parameter_copy = "ALL"

        for index in range(len(prop_all)):
            label = 'COPY - ' + prop_all[index][0]
            
            row = layout.row()
            row.operator("hairmgrplus.copy_parameters", text=label).parameter_copy = str(index)
        

class HAIRMGRPLUS_PT_import_curves(hairmgrplusPanel, bpy.types.Panel):
    bl_label = "Import Curves"
    bl_parent_id = "HAIRMGRPLUS_PT_panel"
    bl_options = {'DEFAULT_CLOSED'}
        
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
        row.operator("hairmgrplus.random_pos", text="Random Pos")          
                    
        row = layout.row()
        row.operator("hairmgrplus.import_from_curves", text="Import from Curves")         
                    
        row = layout.row()        

    
class HAIRMGRPLUS_PT_export_curves(hairmgrplusPanel, bpy.types.Panel):
    bl_label = "Export Curves"
    bl_parent_id = 'HAIRMGRPLUS_PT_panel'
    bl_options = {'DEFAULT_CLOSED'}    
            
    def draw(self, context):
        self.obj = context.object
        layout = self.layout
                
        row = layout.row()        
        row.label(text = "KINDA OFF WORK")
        
        #box = layout.box()            
        row = layout.row()        
        row.prop(getObjLocal(context), 'hamgpExportColl', text = "Export Col.")

        row = layout.row()        
        button = row.operator("hairmgrplus.export", text="Export to Curves")       
        button.export_type = "CURVE"
        
        #row = layout.row()        
        #button = row.operator("hairmgrplus.export", text="Export to Poly") 
        #button.export_type =  "POLY" 

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
    , HAIRMGRPLUS_OT_random_pos
    , HAIRMGRPLUS_OT_import_from_curves
    , HAIRMGRPLUS_OT_export    
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
    