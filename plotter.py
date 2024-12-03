import json
import matplotlib.pyplot as plt

def flatten_2d_array(array_2d):
    """Функция для преобразования двумерного массива в одномерный"""
    return [item for sublist in array_2d for item in sublist]

def plot_graph(filename):
    # Чтение JSON данных из файла
    with open(filename, 'r') as file:
        data = json.load(file)

    # Преобразование двумерных массивов в одномерные
    array_nodes = flatten_2d_array(data['array_nodes'])
    array_f_nodes = flatten_2d_array(data['array_f_nodes'])
    array_x = flatten_2d_array(data['array_x'])
    array_fx = flatten_2d_array(data['array_fx'])
    array_approximation = flatten_2d_array(data['array_approximation'])

    # Построение графика
    plt.figure(figsize =(15,15))
    plt.minorticks_on()
    plt.grid(which='major')
    plt.grid(which='minor')

    # Построение линий
    plt.plot(array_x, array_fx, color='blue', linewidth=2, label='fn')  # Синие линии для исходной функции
    plt.plot(array_x, array_approximation, 'r--', linewidth=1.5, label='L')  # Красные линии для аппроксимации

    # Отображение точек интерполяции
    plt.plot(array_nodes, array_f_nodes, 'gs', markersize=3, label='Points')

    # Настройки осей и легенды
    plt.xlabel('X')
    plt.ylabel('Y')
    plt.legend()
    plt.show()


plot_graph("data.json")