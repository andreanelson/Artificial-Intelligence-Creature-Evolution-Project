# In creature.py
import genome 
from xml.dom.minidom import getDOMImplementation
from enum import Enum
import numpy as np
import math

class MotorType(Enum):
    PULSE = 1
    SINE = 2

class Motor:
    def __init__(self, control_waveform, control_amp, control_freq):
        if control_waveform <= 0.5:
            self.motor_type = MotorType.PULSE
        else:
            self.motor_type = MotorType.SINE
        # self.amp = control_amp
        # self.freq = control_freq
        self.amp = 5
        self.freq = 5
        self.phase = 0

    def get_output(self):
        self.phase = (self.phase + self.freq) % (np.pi * 2)
        if self.motor_type == MotorType.PULSE:
            if self.phase < np.pi:
                output = self.amp
            else:
                output = -self.amp
        if self.motor_type == MotorType.SINE:
            output = np.sin(self.phase)
        return output 
    
class Creature:
    def __init__(self, gene_count):
        self.spec = genome.Genome.get_gene_spec()
        self.dna = genome.Genome.get_random_genome(len(self.spec), gene_count)
        self.flat_links = None
        self.exp_links = None
        self.motors = None
        self.start_position = None
        self.last_position = None
        self.positions = []
        self.previous_distance_to_center = None
        self.fitness = 0

    def get_flat_links(self):
        if self.flat_links == None:
            gdicts = genome.Genome.get_genome_dicts(self.dna, self.spec)
            self.flat_links = genome.Genome.genome_to_links(gdicts)
        return self.flat_links
    
    def get_expanded_links(self):
        self.get_flat_links()
        if self.exp_links is not None:
            return self.exp_links
        
        exp_links = [self.flat_links[0]]
        genome.Genome.expandLinks(self.flat_links[0], 
                                self.flat_links[0].name, 
                                self.flat_links, 
                                exp_links)
        self.exp_links = exp_links
        return self.exp_links

    def to_xml(self):
        self.get_expanded_links()
        domimpl = getDOMImplementation()
        adom = domimpl.createDocument(None, "start", None)
        robot_tag = adom.createElement("robot")
        for link in self.exp_links:
            robot_tag.appendChild(link.to_link_element(adom))
        first = True
        for link in self.exp_links:
            if first:# skip the root node! 
                first = False
                continue
            robot_tag.appendChild(link.to_joint_element(adom))
        robot_tag.setAttribute("name", "pepe") #  choose a name!
        return '<?xml version="1.0"?>' + robot_tag.toprettyxml()

    def get_motors(self):
        self.get_expanded_links()
        if self.motors == None:
            motors = []
            for i in range(1, len(self.exp_links)):
                l = self.exp_links[i]
                m = Motor(l.control_waveform, l.control_amp,  l.control_freq)
                motors.append(m)
            self.motors = motors 
        return self.motors 

    def update_position(self, pos):
        if self.start_position is None:
            self.start_position = pos
        else:
            self.last_position = pos
        self.positions.append(pos)

        # Update fitness based on movement towards or away from the mountain center
        distance_to_center = self.get_distance_to_center()
        if self.previous_distance_to_center is not None:
            if distance_to_center < self.previous_distance_to_center:
                self.fitness += 1  # Reward for moving closer
            else:
                self.fitness -= 1  # Penalty for moving away
        self.previous_distance_to_center = distance_to_center

    def get_distance_travelled(self):
        if self.start_position is None or self.last_position is None:
            return 0
        p1 = np.asarray(self.start_position)
        p2 = np.asarray(self.last_position)
        dist = np.linalg.norm(p1-p2)
        return dist 

    def get_height_reached(self):
        if self.positions:
            heights = [pos[2] for pos in self.positions]
            return max(heights)
        return 0

    def get_stability(self):
        if len(self.positions) < 2:
            return 0
        deviations = [math.sqrt((pos[0] - self.start_position[0])**2 + (pos[1] - self.start_position[1])**2) for pos in self.positions]
        stability = 1 / (1 + np.var(deviations))  # Inverse of variance as a measure of stability
        return stability

    def get_distance_to_center(self):
        if not self.positions:
            return float('inf')  # If no positions, return infinity
        last_position = self.positions[-1]
        center = np.array([0, 0, 0])  # Assuming the mountain center is at the origin
        return np.linalg.norm(np.array(last_position) - center)
        

    # def get_fitness(self):
    #     distance_travelled = self.get_distance_travelled()
    #     height_reached = self.get_height_reached()
    #     stability = self.get_stability()
    #     return distance_travelled + height_reached + stability + self.fitness  # Include the reward/penalty component
    #     # return distance_travelled + height_reached + self.fitness  # Include the reward/penalty component

    def get_fitness(self):
        distance_to_center = self.get_distance_to_center()
        height_reached = self.get_height_reached()
        stability = self.get_stability()
        
        # Inverse distance to center (higher value for closer distance)
        if distance_to_center == 0:
            inverse_distance_to_center = float('inf')  # Avoid division by zero
        else:
            inverse_distance_to_center = 1 / distance_to_center
        
        return stability + height_reached + inverse_distance_to_center + self.fitness  # Include the reward/penalty component



    def reset_position_history(self):
        self.start_position = None
        self.last_position = None
        self.positions = []
        self.previous_distance_to_center = None
        self.fitness = 0

    def update_dna(self, dna):
        self.dna = dna
        self.flat_links = None
        self.exp_links = None
        self.motors = None
        self.start_position = None
        self.last_position = None
        self.positions = []
        self.previous_distance_to_center = None
        self.fitness = 0



# # In creature.py
# import genome 
# from xml.dom.minidom import getDOMImplementation
# from enum import Enum
# import numpy as np
# import math

# class MotorType(Enum):
#     PULSE = 1
#     SINE = 2

# class Motor:
#     def __init__(self, control_waveform, control_amp, control_freq):
#         if control_waveform <= 0.5:
#             self.motor_type = MotorType.PULSE
#         else:
#             self.motor_type = MotorType.SINE
#         self.amp = control_amp
#         self.freq = control_freq
#         self.phase = 0
    

#     def get_output(self):
#         self.phase = (self.phase + self.freq) % (np.pi * 2)
#         if self.motor_type == MotorType.PULSE:
#             if self.phase < np.pi:
#                 output = 1
#             else:
#                 output = -1
            
#         if self.motor_type == MotorType.SINE:
#             output = np.sin(self.phase)
        
#         return output 

# class Creature:
#     def __init__(self, gene_count):
#         self.spec = genome.Genome.get_gene_spec()
#         self.dna = genome.Genome.get_random_genome(len(self.spec), gene_count)
#         self.flat_links = None
#         self.exp_links = None
#         self.motors = None
#         self.start_position = None
#         self.last_position = None
#         self.positions = []

#     def get_flat_links(self):
#         if self.flat_links == None:
#             gdicts = genome.Genome.get_genome_dicts(self.dna, self.spec)
#             self.flat_links = genome.Genome.genome_to_links(gdicts)
#         return self.flat_links
    
#     def get_expanded_links(self):
#         self.get_flat_links()
#         if self.exp_links is not None:
#             return self.exp_links
        
#         exp_links = [self.flat_links[0]]
#         genome.Genome.expandLinks(self.flat_links[0], 
#                                 self.flat_links[0].name, 
#                                 self.flat_links, 
#                                 exp_links)
#         self.exp_links = exp_links
#         return self.exp_links

#     def to_xml(self):
#         self.get_expanded_links()
#         domimpl = getDOMImplementation()
#         adom = domimpl.createDocument(None, "start", None)
#         robot_tag = adom.createElement("robot")
#         for link in self.exp_links:
#             robot_tag.appendChild(link.to_link_element(adom))
#         first = True
#         for link in self.exp_links:
#             if first:# skip the root node! 
#                 first = False
#                 continue
#             robot_tag.appendChild(link.to_joint_element(adom))
#         robot_tag.setAttribute("name", "pepe") #  choose a name!
#         return '<?xml version="1.0"?>' + robot_tag.toprettyxml()

#     def get_motors(self):
#         self.get_expanded_links()
#         if self.motors == None:
#             motors = []
#             for i in range(1, len(self.exp_links)):
#                 l = self.exp_links[i]
#                 m = Motor(l.control_waveform, l.control_amp,  l.control_freq)
#                 motors.append(m)
#             self.motors = motors 
#         return self.motors 
    
#     def update_position(self, pos):
#         if self.start_position == None:
#             self.start_position = pos
#         else:
#             self.last_position = pos
#         self.positions.append(pos)

#     def get_distance_travelled(self):
#         if self.start_position is None or self.last_position is None:
#             return 0
#         p1 = np.asarray(self.start_position)
#         p2 = np.asarray(self.last_position)
#         dist = np.linalg.norm(p1-p2)
#         return dist 

#     def get_height_reached(self):
#         if self.positions:
#             heights = [pos[2] for pos in self.positions]
#             return max(heights)
#         return 0

#     def get_stability(self):
#         if len(self.positions) < 2:
#             return 0
#         deviations = [math.sqrt((pos[0] - self.start_position[0])**2 + (pos[1] - self.start_position[1])**2) for pos in self.positions]
#         stability = 1 / (1 + np.var(deviations))  # Inverse of variance as a measure of stability
#         return stability

#     def get_fitness(self):
#         distance_travelled = self.get_distance_travelled()
#         height_reached = self.get_height_reached()
#         stability = self.get_stability()
#         return distance_travelled + height_reached + stability

#     def reset_position_history(self):
#         self.start_position = None
#         self.last_position = None
#         self.positions = []

#     def update_dna(self, dna):
#         self.dna = dna
#         self.flat_links = None
#         self.exp_links = None
#         self.motors = None
#         self.start_position = None
#         self.last_position = None
#         self.positions = []




# import genome 
# from xml.dom.minidom import getDOMImplementation
# from enum import Enum
# import numpy as np
# import math

# class MotorType(Enum):
#     PULSE = 1
#     SINE = 2

# class Motor:
#     def __init__(self, control_waveform, control_amp, control_freq):
#         if control_waveform <= 0.5:
#             self.motor_type = MotorType.PULSE
#         else:
#             self.motor_type = MotorType.SINE
#         self.amp = control_amp
#         self.freq = control_freq
#         self.phase = 0
    

#     def get_output(self):
#         self.phase = (self.phase + self.freq) % (np.pi * 2)
#         if self.motor_type == MotorType.PULSE:
#             if self.phase < np.pi:
#                 output = 1
#             else:
#                 output = -1
            
#         if self.motor_type == MotorType.SINE:
#             output = np.sin(self.phase)
        
#         return output 

# class Creature:
#     def __init__(self, gene_count):
#         self.spec = genome.Genome.get_gene_spec()
#         self.dna = genome.Genome.get_random_genome(len(self.spec), gene_count)
#         self.flat_links = None
#         self.exp_links = None
#         self.motors = None
#         self.start_position = None
#         self.last_position = None

#     def get_flat_links(self):
#         if self.flat_links == None:
#             gdicts = genome.Genome.get_genome_dicts(self.dna, self.spec)
#             self.flat_links = genome.Genome.genome_to_links(gdicts)
#         return self.flat_links
    
#     def get_expanded_links(self):
#         self.get_flat_links()
#         if self.exp_links is not None:
#             return self.exp_links
        
#         exp_links = [self.flat_links[0]]
#         genome.Genome.expandLinks(self.flat_links[0], 
#                                 self.flat_links[0].name, 
#                                 self.flat_links, 
#                                 exp_links)
#         self.exp_links = exp_links
#         return self.exp_links

#     def to_xml(self):
#         self.get_expanded_links()
#         domimpl = getDOMImplementation()
#         adom = domimpl.createDocument(None, "start", None)
#         robot_tag = adom.createElement("robot")
#         for link in self.exp_links:
#             robot_tag.appendChild(link.to_link_element(adom))
#         first = True
#         for link in self.exp_links:
#             if first:# skip the root node! 
#                 first = False
#                 continue
#             robot_tag.appendChild(link.to_joint_element(adom))
#         robot_tag.setAttribute("name", "pepe") #  choose a name!
#         return '<?xml version="1.0"?>' + robot_tag.toprettyxml()

#     def get_motors(self):
#         self.get_expanded_links()
#         if self.motors == None:
#             motors = []
#             for i in range(1, len(self.exp_links)):
#                 l = self.exp_links[i]
#                 m = Motor(l.control_waveform, l.control_amp,  l.control_freq)
#                 motors.append(m)
#             self.motors = motors 
#         return self.motors 
    
#     def update_position(self, pos):
#         if self.start_position == None:
#             self.start_position = pos
#         else:
#             self.last_position = pos

#     def get_distance_travelled(self):
#         if self.start_position is None or self.last_position is None:
#             return 0
#         p1 = np.asarray(self.start_position)
#         p2 = np.asarray(self.last_position)
#         dist = np.linalg.norm(p1-p2)
#         return dist 

#     def update_dna(self, dna):
#         self.dna = dna
#         self.flat_links = None
#         self.exp_links = None
#         self.motors = None
#         self.start_position = None
#         self.last_position = None