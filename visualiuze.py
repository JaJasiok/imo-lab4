import os
import matplotlib.pyplot as plt

def plot_graph(vertices, file_name):
    # Read edges from file
    with open(os.path.join(os.getcwd(), 'lab4/imo-lab4/' + file_name + '.txt'), 'r') as file:
        lines = file.readlines()

    # Find the index of the empty line separating the two graphs
    empty_line_index = lines.index('\n')

    # Extract edges for the first graph
    edges_graph1 = [tuple(map(int, line.split())) for line in lines[:empty_line_index]]

    # Extract edges for the second graph
    edges_graph2 = [tuple(map(int, line.split())) for line in lines[empty_line_index + 1:] if line.strip()]

    # Plot the vertices
    x = [v[0] for v in vertices]
    y = [v[1] for v in vertices]
    plt.scatter(x, y, color='blue')

    # Plot the edges for the first graph
    for edge in edges_graph1:
        x_values = [vertices[edge[0]][0], vertices[edge[1]][0]]
        y_values = [vertices[edge[0]][1], vertices[edge[1]][1]]
        plt.plot(x_values, y_values, color='red', linewidth=1, alpha=0.5)

    # Plot the edges for the second graph
    for edge in edges_graph2:
        x_values = [vertices[edge[0]][0], vertices[edge[1]][0]]
        y_values = [vertices[edge[0]][1], vertices[edge[1]][1]]
        plt.plot(x_values, y_values, color='green', linewidth=1, alpha=0.5)

    plt.axis('off')  # Turn off the axis
    plt.xticks([])  # Remove x-axis ticks
    plt.yticks([])  # Remove y-axis ticks

    plt.savefig('lab4/imo-lab4/' + file_name + '.png', dpi=100, bbox_inches='tight')
    plt.show()


if __name__ == "__main__":
    # Read vertices from file, skipping first 6 lines
    vertices = []
    with open('lab4/kroB' + '200.tsp', 'r') as file:
        for _ in range(6):
            next(file)  # skip the first 6 lines
        for line in file:
            if line.strip() == 'EOF':  # Stop reading if 'EOF' is encountered
                break
            vertex_id, x, y = map(int, line.split())
            vertices.append((x, y))

    # Get all files with .txt extension in lab4 directory
    txt_files = [file for file in os.listdir('lab4/imo-lab4') if file.endswith('B.txt')]

    for txt_file in txt_files:
        file_name = txt_file.split('.')[0]  # Extract instance name from file name
        plt.clf()
        plot_graph(vertices, file_name)
        plt.clf()