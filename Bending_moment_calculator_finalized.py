import numpy as np
import matplotlib.pyplot as plt
class Shaft:
    def __init__(self, total_length, num_nodes):
        self.total_length = total_length
        self.num_nodes = num_nodes
        # Define a structured data type for nodes (You are NOT creating a TUPLE, this is just the data type definition) Its a fucking array. Then use the node_data to insert - better than console inputs ;)
        self.node_dtype = np.dtype([
            ('position', 'f8'), ('load_z', 'f8'), ('load_y', 'f8'),
            ('load_x', 'f8'), ('dia_gear', 'f8'), ('is_known', 'bool')
        ])
        # Initialize the nodes ARRAY with default values (MF this is an array, dont forget)
        self.nodes = np.zeros(num_nodes, dtype=self.node_dtype)

        self.dictionary_finalized = []

        self.shear_moment = []#Dictionary of shear loads and position. use moment

    def set_node_data(self, node_data):
        if node_data.shape == self.nodes.shape and node_data.dtype == self.nodes.dtype:
            self.nodes = node_data
        else:
            raise ValueError("Node data does not match the expected format or size.")

    def reaction_load_one_unknown(self):
        reaction_load_z = 0
        reaction_load_y = 0
        for node in self.nodes:
            if node['is_known']:
                reaction_load_z += node['load_z']
                reaction_load_y += node['load_y']
        real_reaction_load_z = -reaction_load_z
        real_reaction_load_y = -reaction_load_y

        for node in self.nodes:
            if not node['is_known']:
                node['load_z'] = real_reaction_load_z
                node['load_y'] = real_reaction_load_y
                node['is_known'] = True  # Mark as known after calculation
                break

    def reaction_load_two_unknown(self): 
        sum_of_forces_z = 0
        sum_of_forces_y = 0
        sum_moment_y = 0
        sum_moment_z = 0
        axial_load_moment = int(input('Axial load induces bending about 0: No axial load 1: Z-axis 2: Y-axis: '))

        if axial_load_moment == 0:
            for node in self.nodes:
                if node['is_known']:
                    sum_of_forces_z += node['load_z']
                    sum_of_forces_y += node['load_y']
                sum_moment_y += node['load_y']*node['position']
                sum_moment_z += (-node['load_z']*node['position'])

        elif axial_load_moment == 1:
            axial_moment_direction = int(input('Axial load creates moment about Z-axis. Enter 1(CW) or 2(CCW): '))
            for node in self.nodes:
                if node['is_known']:
                    sum_of_forces_z += node['load_z']
                    sum_of_forces_y += node['load_y']
                if axial_moment_direction == 2:
                    sum_moment_y += node['load_y']*node['position'] + (-node['load_x']*((node['dia_gear'])/2)) #ERROR HERE
                elif axial_moment_direction == 1:
                    sum_moment_y += node['load_y'] * node['position'] + (-node['load_x'] * ((node['dia_gear']) / 2))
                sum_moment_z += (node['load_z']*node['position'])

        else:
            axial_moment_direction = int(input('Axial load creates moment about Y-axis. Enter 1(CW) or 2(CCW): '))
            for node in self.nodes:
                if node['is_known']:
                    sum_of_forces_z += node['load_z']
                    sum_of_forces_y += node['load_y']
                if axial_moment_direction == 2:
                    sum_moment_z += (-node['load_z'] * node['position']) + (-node['load_x'] * ((node['dia_gear']) / 2))
                elif axial_moment_direction == 1:
                    sum_moment_z += (-node['load_z'] * node['position']) + (-node['load_x'] * ((node['dia_gear']) / 2))
                sum_moment_y += node['load_y'] * node['position']

        real_sum_of_forces_z = -sum_of_forces_z
        real_sum_of_forces_y = -sum_of_forces_y
        real_sum_moment_y = -sum_moment_y
        real_sum_moment_z = -sum_moment_z

        unknown_positions = [node['position'] for node in self.nodes if not node['is_known']]
        A = np.array([
            [1, 1],  # Coefficients for forces
            [unknown_positions[0], unknown_positions[1]]  # Coefficients for moments
        ])
        B_z = np.array([real_sum_of_forces_z, real_sum_moment_z])
        B_y = np.array([real_sum_of_forces_y, real_sum_moment_y])

        unknowns_z = np.linalg.solve(A, B_z)
        unknowns_y = np.linalg.solve(A, B_y)

        count = 0
        for node in self.nodes:
            if not node['is_known']:
                node['load_z'] = unknowns_z[count]
                node['load_y'] = unknowns_y[count]
                node['is_known'] = True
                count += 1


#Add all the final values into a blank dictionary then use for diagram. this information is stored ready to be used
    def add_data_to_dictionary(self):
        for node in self.nodes:
            position = node['position']
            load_z = node['load_z']
            load_y = node['load_y']
            load_x = node['load_x']
            dia_gear = node['dia_gear']
            self.add_nodes(position, load_z, load_y, load_x, dia_gear)

    def add_nodes(self, position, load_z, load_y, load_x, dia_gear):
        node_data_finalized = {
            'Position(in)': position,
            'Load z(lb)': load_z,
            'Load y(lb)': load_y,
            'Load x(lb)': load_x,
            'Gear Diameter(in)': dia_gear
        }
        self.dictionary_finalized.append(node_data_finalized)

    def shear_moment_data(self, position, shear_load_z, shear_load_y, bending_moment_z, bending_moment_y):
        shear_moment = {
            'position': position,
            'shear_load_z': shear_load_z,
            'shear_load_y': shear_load_y,
            'bending_moment_z': bending_moment_z,
            'bending_moment_y': bending_moment_y
        }
        self.shear_moment.append(shear_moment)


        def print_shear_moment(self):
            print(f"{'Position':>10} {'Shear Z':>10} {'Shear Y':>10} {'Moment Z':>10} {'Moment Y':>10}")
            for node in self.shear_moment:
                print(
                    f"{node['position']:>10.3f} {node['shear_load_z']:>10.3f} {node['shear_load_y']:>10.3f} {node['bending_moment_z']:>10.3f} {node['bending_moment_y']:>10.3f}")

    def shear_moment_to_list(self):
        shear_load_z = 0
        shear_load_y = 0
        bending_moment_z = 0
        bending_moment_y = 0

        for i, node in enumerate(self.nodes): #This returns indices and each item of the list
            position = node['position']

            if i+1 < len(self.nodes):
                next_node = self.nodes[i+1]
                next_position = next_node['position']
                distance_between_nodes = abs(position - next_position)

                shear_load_z += node['load_z'] #This alows to take one all the shear forces but NOT the last one (they should be 0)
                shear_load_y += node['load_y']

            # Calculate bending moments
                bending_moment_z += shear_load_z * distance_between_nodes
                bending_moment_y += shear_load_y * distance_between_nodes


            # Store the calculated values
            self.shear_moment_data(next_position, shear_load_z, shear_load_y, bending_moment_z, bending_moment_y)

    def shear_diagrams_z(self):
        prev_position = None
        for i, node in enumerate(self.shear_moment):
            if i == 0:
                plt.hlines(y=node['shear_load_z'], xmin=0, xmax=node['position'])
                prev_position = node['position']
            else:
                plt.hlines(y=node['shear_load_z'], xmin=prev_position, xmax=node['position'])
                prev_position = node['position']

            # Annotate each line with its y-value
            plt.annotate(f"{node['shear_load_z']:.3f}",
                         (node['position'] / 2 if i == 0 else (prev_position + node['position']) / 2, node['shear_load_z']), textcoords="offset points", xytext=(0, 10), ha='center')

        plt.title('Shear Force Diagram Z-axis')
        plt.xlabel('Distance (in)')
        plt.ylabel('Shear Force (lb)')
        plt.grid(True)
        plt.show()

    def shear_diagrams_y(self):
        prev_position = None
        for i, node in enumerate(self.shear_moment):
            if i == 0:
                plt.hlines(y=node['shear_load_y'], xmin=0, xmax=node['position'])
                prev_position = node['position']
            else:
                plt.hlines(y=node['shear_load_y'], xmin=prev_position, xmax=node['position'])
                prev_position = node['position']

            # Annotate each line with its y-value
            plt.annotate(f"{node['shear_load_y']:.3f}",(node['position'] / 2 if i == 0 else (prev_position + node['position']) / 2, node['shear_load_y']), textcoords="offset points", xytext=(0, 10), ha='center')

        plt.title('Shear Force Diagram Y-axis')
        plt.xlabel('Distance (in)')
        plt.ylabel('Shear Force (lb)')
        plt.grid(True)
        plt.show()

    def bending_moment_diagram_z(self):
        for node in self.shear_moment:
            plt.plot(node['position'], node['bending_moment_z'], marker='o', linestyle='-', color='red')
            plt.annotate(f"{node['bending_moment_z']:.3f}", (node['position'], node['bending_moment_z']), textcoords="offset points", xytext=(0, 10), ha='center')
        plt.title('Bending Moment Diagram Z-axis')
        plt.xlabel('Distance (in)')
        plt.ylabel('Moment (lb*in)')
        plt.grid(True)
        plt.show()
    def bending_moment_diagram_y(self):
        for node in self.shear_moment:
            plt.plot(node['position'], node['bending_moment_y'], marker='o', linestyle='-', color='red')
            plt.annotate(f"{node['bending_moment_y']:.3f}",(node['position'], node['bending_moment_y']), textcoords="offset points", xytext=(0, 10), ha='center')
        plt.title('Bending Moment Diagram Y-axis')
        plt.xlabel('Distance (in)')
        plt.ylabel('Moment (lb*in)')
        plt.grid(True)
        plt.show()

    def print_shear_moment(self):
        # Header for the table
        print(f"{'Position':>10} {'Shear Z':>10} {'Shear Y':>10} {'Moment Z':>10} {'Moment Y':>10}")

        # Iterate through each node in shear_moment and print the data
        for node in self.shear_moment:
            print(f"{node['position']:>10.3f} "
                  f"{node['shear_load_z']:>10.3f} "
                  f"{node['shear_load_y']:>10.3f} "
                  f"{node['bending_moment_z']:>10.3f} "
                  f"{node['bending_moment_y']:>10.3f}")

#INPUT DATA
total_length = 9.75 #Shaft length
num_nodes = 5      #Total number of nodes
shaft = Shaft(total_length, num_nodes)

# Data Input. Enter ('position, 'load_z', 'load_y', 'load_x', 'dia_gear', 'is_known'). The # of tuples must match num_nodes defined.
node_data = np.array([
    (0, 0, 0, 0, 0, False),  # Node 1, bearing
    (1.25, -24.535, 0, 0, 0, True),   # Node 2 small sprocket
    (8, 0, 0, 0, 0, True), # Node 3: Bigger sprocket
    (9, 3.4076, 12.177971316, 0, 0, True), # Node 4: bearing
    (9.75, 0, -11.240447155, 0, 0, True),# Node 5: Spline
], dtype=shaft.node_dtype) #Data type assignment


shaft.set_node_data(node_data) #Set self.nodes with input data
count_unknown = 0
for unknown in shaft.nodes:
    if unknown['is_known'] == False:
        count_unknown += 1
if count_unknown == 1:
    shaft.reaction_load_one_unknown()
elif count_unknown == 2:
    shaft.reaction_load_two_unknown()


shaft.add_data_to_dictionary()

print("\n")

print("Node Data after Reaction Load Calculation:")
for node in shaft.nodes:
    print(f"Position: {node['position']:.3f}, Load Z: {node['load_z']:.3f}, Load Y: {node['load_y']:.3f}, Load X: {node['load_x']:.3f}, Gear Diameter: {node['dia_gear']:.3f}, Known: {node['is_known']}")

shaft.shear_moment_to_list()
print("\n")
print(f"{'Position':>10} {'Shear Z':>10} {'Shear Y':>10} {'Moment Z':>10} {'Moment Y':>10}")
for i, node in enumerate(shaft.shear_moment):
    # Check if it's not the last node or if it is different from the next node
    if i < len(shaft.shear_moment) - 1 and shaft.shear_moment[i] != shaft.shear_moment[i + 1]:
        print(f"{node['position']:>10.3f} "
              f"{node['shear_load_z']:>10.3f} "
              f"{node['shear_load_y']:>10.3f} "
              f"{node['bending_moment_z']:>10.3f} "
              f"{node['bending_moment_y']:>10.3f}")
    elif i == len(shaft.shear_moment) - 1:  # Always print the last node
        print(f"{node['position']:>10.3f} "
              f"{node['shear_load_z']:>10.3f} "
              f"{node['shear_load_y']:>10.3f} "
              f"{node['bending_moment_z']:>10.3f} "
              f"{node['bending_moment_y']:>10.3f}")



shaft.shear_diagrams_z()
shaft.shear_diagrams_y()
shaft.bending_moment_diagram_z()
shaft.bending_moment_diagram_y()

#Takeaway
#MODULARITY MF
#USE fucking datatypes (numpy)
#Use FAT fucking dictionaries, no small mfs. With enumerate this shit is a handbook
#Name better
#VALIDATION AT EACH DEF OR U RESTART
#MOST IMPORTANT: STEP BY STEP PROCESS. USE FUCKIONG ENUMERATE (INDICES and LIST ITERATION at the same time)
