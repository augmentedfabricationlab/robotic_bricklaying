from __future__ import print_function
from __future__ import absolute_import
from __future__ import division

from compas.datastructures import Datastructure
from compas.datastructures import Graph
from compas.datastructures import AssemblyError
from compas.geometry import Frame, Translation, Rotation, Transformation, Point, Vector, Plane
from compas_rhino.conversions import plane_to_compas_frame, point_to_compas

from assembly_information_model import Assembly
from .part import CAEPart as Part

import math
import Rhino.Geometry as rg


class CAEAssembly(Assembly):
    """A data structure for managing the connections between different parts of an assembly.

    Parameters
    ----------
    name : str, optional
        The name of the assembly.
    **kwargs : dict, optional
        Additional keyword arguments, which are stored in the attributes dict.

    Attributes
    ----------
    graph : :class:`compas.datastructures.Graph`
        The graph that is used under the hood to store the parts and their connections.

    See Also
    --------
    :class:`compas.datastructures.Graph`
    :class:`compas.datastructures.Mesh`
    :class:`compas.datastructures.VolMesh`

    """

    def __init__(self, name=None, **kwargs):
        super(CAEAssembly, self).__init__()

    def export_to_json(self, path, is_built=False):

        # TODO!
        self.graph.update_default_node_attributes({"is_built":False})
        for key in self.parts():
            self.graph.node_attribute(key, "is_built", is_built)

        self.to_json(path)



    def add_to_assembly(self,
                        brick_type, 
                        brick_geometry, 
                        brick_insulated,
                        fixed = True, 
                        frame=None
                        ): 
        """Create a brick with a specified type and to add it to the assembly"""

        if frame is None:
            frame = frame

        if brick_type == "full":
            brick = brick_geometry

        if brick_type == "insulated":
            brick = brick_insulated
        transformed_brick = brick.transformed(Transformation.from_frame(frame))
        brick_height = brick_geometry.shape.zsize
        gripping_frame = frame.transformed(Translation.from_vector(frame.zaxis*((brick_height-0.020)/2)))
        R = Rotation.from_axis_and_angle(gripping_frame.xaxis, math.radians(180), gripping_frame.point)
        gripping_frame.transform(R)
        transformed_brick.gripping_frame = gripping_frame

        self.add_part(transformed_brick, attr_dict={"brick_type": brick_type, "fixed": fixed})

    def generate_french_bond(self,
                    brick_geometry,
                    brick_insulated,
                    initial_brick_center,
                    line_length,          
                    plane,
                    course_is_odd,
                    j):

        brick_spacing = 0.015
        brick_width = brick_geometry.shape.xsize
        brick_length = brick_geometry.shape.ysize

        params = {
            "brick_geometry": brick_geometry, 
            "brick_insulated": brick_insulated
        }
        
        center_frame = plane_to_compas_frame(plane)
        num_bricks1 = math.floor(line_length / (((brick_width+2*brick_length) + 3*brick_spacing)))

        for i in range(num_bricks1):
            T = plane.XAxis * -(i*(2*(brick_spacing+brick_length)))
            translation = Translation.from_vector(T)
            
            # Apply translation to the initial brick center
            brick_center = initial_brick_center + T
            brick_frame = Frame(point_to_compas(brick_center), center_frame.xaxis, center_frame.yaxis)
            
            # Transform the frame with translation
            current_frame = brick_frame.transformed(translation)
            # Add the brick to the assembly
            if course_is_odd:
                self.add_to_assembly(brick_type="full", fixed=False, frame=current_frame, **params)
            else:
                T = plane.XAxis * -((((brick_length+brick_spacing+brick_width)/2)))
                translation = Translation.from_vector(T)
                current_frame = current_frame.transformed(translation)
                self.add_to_assembly(brick_type="full", fixed=False, frame=current_frame, **params)

            
            T = plane.XAxis * -((((brick_length+brick_spacing)*1.5)))
            T2 = plane.XAxis * -(i*(((2*(brick_length+brick_spacing)))))
            T1 = plane.YAxis * ((brick_width-brick_length)/2)
            
            translation = Translation.from_vector(T)
            Translation2 = Translation.from_vector(T1)
            Translation3 = Translation.from_vector(T2)
            
            # Create a rotation transformation (90 degrees around Z-axis)
            R = Rotation.from_axis_and_angle(current_frame.zaxis, math.radians(90), brick_frame.point)
            
            # Apply rotation
            rotated_frame = brick_frame.transformed(R)
            
            # Apply translation to the rotated frame
            current_frame = rotated_frame.transformed(translation*Translation2*Translation3)
            
            # Add the rotated brick to the assembly
            if course_is_odd:
                self.add_to_assembly(brick_type="full", fixed=True, frame=current_frame, **params)
            else:
                T = plane.XAxis * -((((brick_length+brick_spacing+brick_width)/2)))
                translation = Translation.from_vector(T)
                current_frame = current_frame.transformed(translation)
                self.add_to_assembly(brick_type="full", fixed=True, frame=current_frame, **params)

            T = plane.XAxis * -((((brick_length+brick_spacing))))
            translation = Translation.from_vector(T)

            current_frame = current_frame.transformed(translation)

            self.add_to_assembly(brick_type="full", fixed=True, frame=current_frame, **params)





    def generate_cross_bond(self,
                            brick_geometry,
                            brick_insulated,
                            initial_brick_center,
                            line_length,          
                            plane,
                            course_is_odd,
                            j):
        """
        Generates a Cross Bond pattern of bricks.

        Parameters:
        - brick_geometry: Geometry of the brick.
        - brick_insulated: Insulated brick geometry.
        - initial_brick_center: Starting center point for brick placement (COMPAS Point).
        - num_bricks: Number of bricks in the course.
        - plane: Reference plane for brick placement.
        - course_is_odd: Boolean indicating if the course is odd.
        """
        brick_spacing = 0.015
        brick_width = brick_geometry.shape.xsize
        brick_length = brick_geometry.shape.ysize

        params = {
            "brick_geometry": brick_geometry, 
            "brick_insulated": brick_insulated
        }
        
        center_frame = plane_to_compas_frame(plane)
        ornament = 0

        if course_is_odd:
            num_bricks1 = math.floor(line_length / (brick_width+brick_spacing))
            num_bricks2 = math.floor(line_length / (brick_length+brick_spacing))
            # Odd courses: Bricks laid with the long side facing out
            for i in range(num_bricks1):
                # Calculate translation vector for the current brick
                T = plane.XAxis * -(i * (((brick_width+ brick_spacing)/2)))
                translation = Translation.from_vector(T)
                
                # Apply translation to the initial brick center
                brick_center = initial_brick_center + T
                brick_frame = Frame(point_to_compas(brick_center), center_frame.xaxis, center_frame.yaxis)
                
                # Transform the frame with translation
                current_frame = brick_frame.transformed(translation)
                # Add the brick to the assembly
                if ornament == 0 or ornament ==1:
                    self.add_to_assembly(brick_type="full", fixed=True, frame=current_frame, **params)

                T1 = plane.YAxis * ((brick_width+ brick_length + (2*brick_spacing)))
                T2 = plane.XAxis * -((brick_length + brick_spacing)/2)
                translation1 = Translation.from_vector(T1)
                translation2 = Translation.from_vector(T2)
                current_frame = current_frame.transformed(translation1*translation2)

    
                self.add_to_assembly(brick_type="insulated", fixed=True, frame=current_frame, **params)







            for i in range(num_bricks2):
                T = plane.XAxis * -(i * (((brick_length+brick_spacing)/2)))
                T1 = plane.YAxis * (((brick_width-brick_length)/2)+ (brick_length+brick_spacing))
                
                translation = Translation.from_vector(T)
                Translation2 = Translation.from_vector(T1)
                
                # Create the initial brick frame
                brick_center = initial_brick_center + T
                brick_frame = Frame(point_to_compas(brick_center), center_frame.xaxis, center_frame.yaxis)
                
                # Create a rotation transformation (90 degrees around Z-axis)
                R = Rotation.from_axis_and_angle(brick_frame.zaxis, math.radians(90), brick_frame.point)
                
                # Apply rotation
                rotated_frame = brick_frame.transformed(R)
                
                # Apply translation to the rotated frame
                current_frame = rotated_frame.transformed(translation*Translation2)
                
    

                self.add_to_assembly(brick_type="full", fixed=True, frame=current_frame, **params)



                        
        elif not course_is_odd:
            num_bricks = math.floor(line_length / (brick_length+brick_spacing))
            # Even courses: Bricks laid with the short side facing out (rotated by 90 degrees)
            for i in range(num_bricks):
                # Calculate translation vector based on brick length
                T = plane.XAxis * -(i * (((brick_length+brick_spacing)/2)))
                T1 = plane.YAxis * ((brick_width-brick_length)/2)
                
                translation = Translation.from_vector(T)
                Translation2 = Translation.from_vector(T1)
                
                # Create the initial brick frame
                brick_center = initial_brick_center + T
                brick_frame = Frame(point_to_compas(brick_center), center_frame.xaxis, center_frame.yaxis)
                
                # Create a rotation transformation (90 degrees around Z-axis)
                R = Rotation.from_axis_and_angle(brick_frame.zaxis, math.radians(90), brick_frame.point)
                
                # Apply rotation
                rotated_frame = brick_frame.transformed(R)
                
                # Apply translation to the rotated frame
                current_frame = rotated_frame.transformed(translation*Translation2)
                
                # Add the rotated brick to the assembly
  

                if ornament == 0: 
                    if i % 2 == 0 and j% 4 == 0:
                        self.add_to_assembly(brick_type="full", fixed=False, frame=current_frame, **params)
                    elif i % 2 != 0 and j% 4 != 0:
                        self.add_to_assembly(brick_type="full", fixed=False, frame=current_frame, **params)
                    else:
                        self.add_to_assembly(brick_type="full", fixed=True, frame=current_frame, **params)

                if ornament == 1: 
                    if i % 2 == 0:
                        self.add_to_assembly(brick_type="full", fixed=False, frame=current_frame, **params)
                    else:
                        self.add_to_assembly(brick_type="full", fixed=True, frame=current_frame, **params)


                if ornament == 3: 
                    if i % 2 == 0:
                        self.add_to_assembly(brick_type="full", fixed=False, frame=current_frame, **params)
                    else:
                        self.add_to_assembly(brick_type="full", fixed=True, frame=current_frame, **params)


                T2 = plane.YAxis * (brick_width+ brick_spacing)
                T3 = plane.XAxis * ((brick_length+brick_spacing)/2)
                Translation3 = Translation.from_vector(T2)
                Translation4= Translation.from_vector(T3)
                current_frame = current_frame.transformed(Translation3 * Translation4)
                self.add_to_assembly(brick_type="insulated", fixed=True, frame=current_frame, **params)
  

    def generate_flemish_bond(self,
                                brick_geometry,
                                brick_insulated,
                                initial_brick_center,
                                num_bricks,
                                plane,
                                course_is_odd,
                                ):

        brick_spacing = 0.015
        #brick_height = brick_geometry.shape.zsize
        brick_width = brick_geometry.shape.xsize
        brick_length = brick_geometry.shape.ysize
        #course = brick_height+brick_spacing

        params = {"brick_geometry": brick_geometry, 
                    "brick_insulated": brick_insulated
                    }

        center_frame = plane_to_compas_frame(plane)

        for i in range(num_bricks):


            T = plane.XAxis * (i * ((brick_width +brick_length)/2 + brick_spacing))
            #Shifting every second row
            if course_is_odd == True:
                T += plane.XAxis * (brick_length + brick_width + brick_spacing)/2
            brick_center = initial_brick_center + T
            brick_frame = Frame(point_to_compas(brick_center), center_frame.xaxis, center_frame.yaxis)
            
            if i % 2 == 0:

                #self-shading bricks - first layer
                R = Rotation.from_axis_and_angle(brick_frame.zaxis, math.radians(90), point=brick_frame.point)  # Rotation
                T1 = Translation.from_vector(brick_frame.xaxis * (brick_length + brick_spacing)/2)
                current_frame = brick_frame.transformed(R*T1)
                self.add_to_assembly(brick_type = "full", fixed = False, frame = current_frame, **params) 

                #bricks - second layer - insulated
                T2 = Translation.from_vector(brick_frame.yaxis * (brick_width + brick_spacing))
                current_frame = brick_frame.transformed(T2)
                self.add_to_assembly(brick_type = "insulated", fixed = True, frame = current_frame, **params) 

                #bricks - second layer - insulated
                T4 = Translation.from_vector(current_frame.yaxis * (brick_length + brick_spacing))
                current_frame = current_frame.transformed(T4)
                self.add_to_assembly(brick_type = "insulated", fixed = True, frame = current_frame, **params) 

            else:
                #bricks - first layer - solid
                self.add_to_assembly(brick_type = "full", fixed = True, frame = brick_frame, **params) 

                #bricks - first layer - insulated
                T1 = Translation.from_vector(brick_frame.yaxis * (brick_length + brick_spacing))
                current_frame = brick_frame.transformed(T1)
                self.add_to_assembly(brick_type = "insulated", fixed = True, frame = current_frame , **params) 

                #bricks - second layer
                T2 = Translation.from_vector(brick_frame.yaxis * ((brick_length+brick_width)/2 + brick_spacing))
                R = Rotation.from_axis_and_angle(current_frame.zaxis, math.radians(90), point=current_frame.point)
                current_frame = current_frame.transformed(T2*R)
                self.add_to_assembly(brick_type = "insulated", fixed = True, frame = current_frame, **params)

            #Outter frame of the wall:

            if course_is_odd:

                if i == 0:
                    T = Translation.from_vector(-brick_frame.xaxis * (brick_length/2 + brick_width/2 +brick_spacing))
                    current_frame = brick_frame.transformed(T)
                    self.add_to_assembly(brick_type = "full", fixed = True, frame = current_frame, **params) 

                    T = Translation.from_vector(current_frame.yaxis * (brick_length + brick_spacing))
                    current_frame = current_frame.transformed(T)
                    self.add_to_assembly(brick_type = "insulated", fixed = True, frame = current_frame, **params) 

                    T = Translation.from_vector(brick_frame.xaxis * (brick_width/2 + brick_length/2 + brick_spacing))
                    R = Rotation.from_axis_and_angle(current_frame.zaxis, math.radians(90), point=current_frame.point) 
                    current_frame = current_frame.transformed(R*T)
                    self.add_to_assembly(brick_type = "insulated", fixed = True, frame = current_frame, **params) 


    def generate_other_bond(self,
                    brick_geometry,
                    brick_insulated,
                    initial_brick_center,
                    num_bricks,
                    plane,
                    course_is_odd,
                    ):

        brick_spacing = 0.015
        #brick_height_f = brick_geometry.shape.zsize
        #brick_width_f = brick_geometry.shape.xsize
        brick_length_f = brick_geometry.shape.ysize
        brick_width_i = brick_insulated.shape.xsize
        brick_length_i = brick_insulated.shape.ysize
      
        

        params = {"brick_geometry": brick_geometry, 
                    "brick_insulated": brick_insulated 
                    }

        center_frame = plane_to_compas_frame(plane)

        for i in range(num_bricks):
            T = plane.XAxis * (i * (((5*brick_length_i)/3)+0.002))
            #Shifting every second row
            if course_is_odd == True:
                T += plane.XAxis * ((brick_length_i+brick_width_i+ 3* brick_spacing)/2-0.003)
            brick_center = initial_brick_center + T

            brick_frame = Frame(point_to_compas(brick_center), center_frame.xaxis, center_frame.yaxis)

            if not course_is_odd:

                if i % 2 != 0:
                    

                    #self-shading bricks - first layer
                    R = Rotation.from_axis_and_angle(brick_frame.zaxis, math.radians(90), point=brick_frame.point)  
                    T1 = Translation.from_vector(brick_frame.xaxis * (brick_length_f + brick_spacing)/2)
                    current_frame = brick_frame.transformed(R*T1)

                    self.add_to_assembly(brick_type = "full", fixed = False, frame = current_frame, **params) 

                    T2 = Translation.from_vector(current_frame.xaxis * (brick_width_i + brick_spacing))
                    current_frame2 = current_frame.transformed(T2)
                    self.add_to_assembly(brick_type = "insulated", fixed = True, frame = current_frame2, **params) 

                else:
                    #bricks - first layer - solid
                    self.add_to_assembly(brick_type = "full", fixed = True, frame = brick_frame, **params) 

                    T1 = Translation.from_vector(brick_frame.yaxis * (brick_length_i + brick_width_i + 2 * brick_spacing))
                    current_frame = brick_frame.transformed(T1)
                    T12 = Translation.from_vector(current_frame.yaxis * (-0.003))
                    current_frame3 = current_frame.transformed(T12)
                    self.add_to_assembly(brick_type = "insulated", fixed = True, frame = current_frame3 , **params) 

                    R = Rotation.from_axis_and_angle(current_frame.zaxis, math.radians(90), point=current_frame.point)
                    T1 = Translation.from_vector(current_frame.yaxis * (((brick_length_i + (brick_spacing))/2)))
                    current_frame = current_frame.transformed(R*T1)
                    T2 = Translation.from_vector(current_frame.xaxis * -((brick_width_i + brick_length_i + (2.5 * brick_spacing+0.002))/2))
                    current_frame = current_frame.transformed(T2)
                    self.add_to_assembly(brick_type = "insulated", fixed = True, frame = current_frame , **params)

                    T3 = Translation.from_vector(current_frame.yaxis * -(brick_width_i + (2*brick_spacing-0.006))/2)
                    current_frame = current_frame.transformed(T3)
                    self.add_to_assembly( brick_type = "insulated", fixed = True, frame = current_frame , **params)



            #Outter frame of the wall:
        
            if course_is_odd:

                if i == 0:
                    #self.add_to_assembly(brick_type = "full", fixed = True, frame = brick_frame, **params) 
                    R = Rotation.from_axis_and_angle(brick_frame.zaxis, math.radians(90), point=brick_frame.point)  
                    current_frame = brick_frame.transformed(R)
                    T = Translation.from_vector(current_frame.xaxis *((brick_length_f + brick_spacing)/2))
                    current_frame = current_frame.transformed(T)
                    T2 = Translation.from_vector(current_frame.yaxis * (brick_length_i + brick_length_i/2 + brick_width_i + 3 * brick_spacing)/2)
                    current_frame = current_frame.transformed(T2)
                    self.add_to_assembly(brick_type = "full", fixed = False, frame = current_frame, **params) 

                    T = Translation.from_vector(current_frame.xaxis *((brick_width_i+brick_spacing)))
                    current_frame = current_frame.transformed(T)
                    self.add_to_assembly(brick_type = "insulated", fixed = True, frame = current_frame, **params) 

                if i >= 0 and i < num_bricks-1:


                    if i % 2 != 0:
                        #self.add_to_assembly(brick_type = "full", fixed = True, frame = current_frame, **params) 
                
                        #self-shading bricks - first layer
                        R = Rotation.from_axis_and_angle(brick_frame.zaxis, math.radians(90), point=brick_frame.point)  
                        T1 = Translation.from_vector(brick_frame.xaxis * (brick_length_f + brick_spacing)/2)
                        current_frame = brick_frame.transformed(R*T1)

                        self.add_to_assembly(brick_type = "full", fixed = False, frame = current_frame, **params) 

                        T2 = Translation.from_vector(current_frame.xaxis * (brick_width_i + brick_spacing))
                        current_frame2 = current_frame.transformed(T2)
                        self.add_to_assembly(brick_type = "insulated", fixed = True, frame = current_frame2, **params) 

                    else:
                        #bricks - first layer - solid
                        self.add_to_assembly(brick_type = "full", fixed = True, frame = brick_frame, **params) 

                        T1 = Translation.from_vector(brick_frame.yaxis * (brick_length_i + brick_width_i + 2 * brick_spacing))
                        current_frame = brick_frame.transformed(T1)
                        T12 = Translation.from_vector(current_frame.yaxis * (-0.003))
                        current_frame3 = current_frame.transformed(T12)
                        self.add_to_assembly(brick_type = "insulated", fixed = True, frame = current_frame3 , **params) 

                        R = Rotation.from_axis_and_angle(current_frame.zaxis, math.radians(90), point=current_frame.point)
                        T1 = Translation.from_vector(current_frame.yaxis * (((brick_length_i + (brick_spacing))/2)))
                        current_frame = current_frame.transformed(R*T1)
                        T2 = Translation.from_vector(current_frame.xaxis * -((brick_width_i + brick_length_i + (2.5 * brick_spacing+0.002))/2))
                        current_frame = current_frame.transformed(T2)
                        self.add_to_assembly(brick_type = "insulated", fixed = True, frame = current_frame , **params)

                        T3 = Translation.from_vector(current_frame.yaxis * -(brick_width_i + (2*brick_spacing-0.006))/2)
                        current_frame = current_frame.transformed(T3)
                        self.add_to_assembly( brick_type = "insulated", fixed = True, frame = current_frame , **params)


                elif i == num_bricks-1:

                    R = Rotation.from_axis_and_angle(brick_frame.zaxis, math.radians(90), point=brick_frame.point)  
                    current_frame = brick_frame.transformed(R)
                    T = Translation.from_vector(current_frame.xaxis *((brick_length_f + brick_spacing)/2))
                    current_frame = current_frame.transformed(T)
                    T2 = Translation.from_vector(current_frame.yaxis * - (brick_length_i/2)/2)
                    current_frame = current_frame.transformed(T2)
                    self.add_to_assembly(brick_type = "full", fixed = False, frame = current_frame, **params) 

                    T = Translation.from_vector(current_frame.xaxis * ((brick_width_i+brick_spacing)))
                    current_frame = current_frame.transformed(T)
                    self.add_to_assembly(brick_type = "insulated", fixed = True, frame = current_frame, **params)
             




    def generate_wall(self,
                            bond,
                            brick_geometry,
                            brick_insulated,
                            plane,
                            lines):
        brick_spacing = 0.015
        brick_height = brick_geometry.shape.zsize
        brick_width = brick_geometry.shape.xsize
        brick_length = brick_geometry.shape.ysize

        center_frame = plane_to_compas_frame(plane)
        total_length = 0
        start_x = lines[1].FromX

        params = {"brick_geometry": brick_geometry,
                  "brick_insulated": brick_insulated,
                  "plane": plane
        }

        for j, line in enumerate(lines):
            line_length = line.Length
            halfway_point = line.PointAt(line.Length / 2)
            num_bricks = math.floor(line_length / ((brick_width+brick_length)/2+brick_spacing))

            if num_bricks % 2 == 0:
                num_bricks = num_bricks-1

            course_is_odd = j %2 != 0 #check if course is odd or even

            if course_is_odd == True and num_bricks % 2 != 0:
                num_bricks = num_bricks-1

            # Making sure all the lines are the same
            start_point = line.From
            start_point.X = start_x 
            line.From = start_point
            initial_brick_center = line.From


            T = plane.XAxis * (brick_length+brick_spacing)
            initial_brick_center += T
        

            if j == 0:
                #Calculating the length of the wall
                half_bricks = math.ceil(num_bricks / 2)
                total_length += half_bricks * (brick_length+brick_spacing) + (num_bricks-half_bricks)*(brick_width+brick_spacing)
                total_length += 2 * (brick_length/2) 

            #Pick the bond   
            if bond == 0:
                self.generate_flemish_bond(
                                    initial_brick_center = initial_brick_center,
                                    num_bricks = num_bricks,
                                    course_is_odd = course_is_odd,
                                    **params) 
            if bond == 1:
                self.generate_other_bond(
                        initial_brick_center = initial_brick_center,
                        num_bricks = num_bricks,
                        course_is_odd = course_is_odd,
                        **params) 
            if bond == 2:
                self.generate_other_bond_archive(
                        initial_brick_center = initial_brick_center,
                        num_bricks = num_bricks,
                        course_is_odd = course_is_odd,
                        **params) 
                
            if bond == 3:
                self.generate_cross_bond(
                        initial_brick_center = initial_brick_center,
                        line_length= line_length,
                        course_is_odd = course_is_odd,
                        j=j,
                        **params) 

            if bond == 4:
                self.generate_french_bond(
                        initial_brick_center = initial_brick_center,
                        line_length= line_length,
                        course_is_odd = course_is_odd,
                        j=j,
                        **params)
                
                
        return total_length
    
    def apply_gradient(self, transform, values, keys):
        i = 0
        for key in keys:
            print(key)
            if key == 4:
                pass
            else:
                
                part = self.part(key)

                y_translation = values[i]
                y_translation *= -0.08

                
                translation_vector = part.frame.xaxis * y_translation
    
                    
                translation = Translation.from_vector(translation_vector)

                # Apply the translation to the frame
                transformed_frame = part.frame.transformed(translation)

                # Update the geometry position
                part.transform(translation)
            i = i +1


    def generate_other_bond_archive(self,
                    brick_geometry,
                    brick_insulated,
                    initial_brick_center,
                    num_bricks,
                    plane,
                    course_is_odd,
                    ):

        brick_spacing = 0.015
        #brick_height_f = brick_geometry.shape.zsize
        brick_width_f = brick_geometry.shape.xsize
        brick_length_f = brick_geometry.shape.ysize

        brick_width_i = brick_insulated.shape.xsize
        brick_length_i = brick_insulated.shape.ysize
      
        

        params = {"brick_geometry": brick_geometry, 
                    "brick_insulated": brick_insulated 
                    }

        center_frame = plane_to_compas_frame(plane)

        for i in range(num_bricks):
            T = plane.XAxis * (i * ((brick_width_f +brick_length_f)/2 + brick_spacing))
            #Shifting every second row
            if course_is_odd == True:
                T += plane.XAxis * (brick_length_f + brick_width_f + brick_spacing)/2
            brick_center = initial_brick_center + T

            brick_frame = Frame(point_to_compas(brick_center), center_frame.xaxis, center_frame.yaxis)

            if i % 2 != 0:

                #self-shading bricks - first layer
                R = Rotation.from_axis_and_angle(brick_frame.zaxis, math.radians(90), point=brick_frame.point)  
                T1 = Translation.from_vector(brick_frame.xaxis * (brick_length_f + brick_spacing)/2)
                current_frame = brick_frame.transformed(R*T1)

                self.add_to_assembly(brick_type = "full", fixed = False, frame = current_frame, **params) 

                T2 = Translation.from_vector(current_frame.xaxis * (brick_width_f + brick_spacing))
                current_frame2 = current_frame.transformed(T2)
                self.add_to_assembly(brick_type = "insulated", fixed = True, frame = current_frame2, **params) 

            else:
                #bricks - first layer - solid
                #T10 = Translation.from_vector(brick_frame.xaxis * (0.035))
                #current_frame = brick_frame.transformed(T10)
                self.add_to_assembly(brick_type = "full", fixed = True, frame = brick_frame, **params) 
       


                T1 = Translation.from_vector(brick_frame.yaxis * (brick_length_f + brick_width_f+ 2 * brick_spacing))
                current_frame = brick_frame.transformed(T1)
                self.add_to_assembly(brick_type = "insulated", fixed = True, frame = current_frame , **params) 

                R = Rotation.from_axis_and_angle(current_frame.zaxis, math.radians(90), point=current_frame.point)
                T1 = Translation.from_vector(current_frame.yaxis * (((brick_length_f + brick_spacing)/2)))
                current_frame = current_frame.transformed(R*T1)
                T2 = Translation.from_vector(current_frame.xaxis * -((brick_width_f + brick_length_f + 2*brick_spacing)/2))
                current_frame = current_frame.transformed(T2)
                self.add_to_assembly(brick_type = "insulated", fixed = True, frame = current_frame , **params)

                T3 = Translation.from_vector(current_frame.yaxis * -(brick_width_f + brick_spacing)/2)
                current_frame = current_frame.transformed(T3)
                self.add_to_assembly( brick_type = "insulated", fixed = True, frame = current_frame , **params)



            #Outter frame of the wall:
        
            if course_is_odd:

                if i == 0:
                    R = Rotation.from_axis_and_angle(brick_frame.zaxis, math.radians(90), point=brick_frame.point)  
                    current_frame = brick_frame.transformed(R)
                    T = Translation.from_vector(current_frame.xaxis *((brick_length_f+brick_spacing)/2))
                    current_frame = current_frame.transformed(T)
                    T2 = Translation.from_vector(current_frame.yaxis * (brick_length_f + brick_width_f + 2* brick_spacing)/2)
                    current_frame = current_frame.transformed(T2)
                    self.add_to_assembly(brick_type = "full", fixed = False, frame = current_frame, **params) 

                    T = Translation.from_vector(current_frame.xaxis *((brick_width_f+brick_spacing)))
                    current_frame = current_frame.transformed(T)
                    self.add_to_assembly(brick_type = "insulated", fixed = True, frame = current_frame, **params) 





