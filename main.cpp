#include <iostream>
#include <cstdlib>
#include "best_approximation.hpp"
using namespace std;

int main(){
    // Создание основных массивов
    auto* array_global_nodes = gen_uniform_grid(Task_const::H, Task_const::M); // Генерируем глобальную сетку
    auto** array_2d_nodes = gen_2d_arr_uniform(array_global_nodes); // Разбиваем глобальную сетку на подсетки (с перекрытием в 1 узел)
    auto** array_2d_random = gen_random_2d_arr_in_local(array_2d_nodes); // Генерируем массив случайных точек на каждом КЭ, (В нем содержатся только случаные точки, без узлов)
    auto** array_2d_func_in_rand = gen_func_2d_arr(array_2d_random, Task_const::L); // Применяем заданную функцию к массиву случаных точек
    auto** array_2d_func_in_nodes = gen_func_2d_arr(array_2d_nodes,Task_const::N);
    
    // Создание и решение слау
    auto** global_matrix = gen_global_matrix(array_2d_nodes, array_2d_random); // Собираем глобальную матрицу для решения слау
    auto* global_b_vector = gen_global_vector_b(array_2d_nodes, array_2d_random, array_2d_func_in_rand); // Собираем глобальную правую часть для решения слау
    auto** coefficients = solve_system(global_matrix, global_b_vector); // Решаем слау, и находим коэффициенты для аппрокцисмации
    
    // Готовим массив rand для построения на нем наилучшего приближения
    auto** arr_approx_in_rand_points = best_approximation_global(coefficients, array_2d_nodes, array_2d_random, Task_const::L); // На основе найденных коэффициентов и базисных функций строим аппроксимацию
    
    // Готовим сетку h/100
    std::size_t count_h_100_points_in_ever_elem; //Количество точек с шагом h/100 на каждом элементе
    auto** array_2d_h_100_points = gen_2d_uniform_between_nodes(false, array_2d_nodes, count_h_100_points_in_ever_elem, Task_const::STEP_H_100); // Создаем сетку с шагом h_100
    
    auto** arr_approx_in_h_100_points = best_approximation_global(coefficients, array_2d_nodes, array_2d_h_100_points, count_h_100_points_in_ever_elem); // Вычисляем наилучшее приближение в сетке h/100
    // Вычисляем значения функции в сетке h_100
    
    auto** array_func_in_h_100_points = gen_func_2d_arr(array_2d_h_100_points, count_h_100_points_in_ever_elem);
    //Вычисление погрешностей
    auto errors_random_points = calculate_errors(array_2d_func_in_rand, arr_approx_in_rand_points, Task_const::L, Task_const::K); // Вычисляем погрешности в случайных точках
    auto errors_h_100_points = calculate_errors(array_func_in_h_100_points, arr_approx_in_h_100_points, count_h_100_points_in_ever_elem, Task_const::K); // Вычисляем погрешности в точках h/100

    // Готовим сетку M_viz
    std::size_t count_M_viz_points_in_ever_elem;
    auto** array_M_viz_points = gen_2d_uniform_between_nodes(true, array_2d_nodes, count_M_viz_points_in_ever_elem, Task_const::STEP_M_viz); // Уплотняем существующий массив с сеткой для построения графика
    auto** arr_approx_in_M_viz_points = best_approximation_global(coefficients, array_2d_nodes, array_M_viz_points, count_M_viz_points_in_ever_elem);
    auto** arr_2d_func_in_M_viz_points = gen_func_2d_arr(array_M_viz_points, count_M_viz_points_in_ever_elem);
    
    write_data_to_file(array_2d_nodes, array_2d_func_in_nodes, array_M_viz_points, arr_2d_func_in_M_viz_points, arr_approx_in_M_viz_points, count_M_viz_points_in_ever_elem);
    print_error_table(errors_random_points, errors_h_100_points);
    system("python plotter.py");

    delete[] array_global_nodes;
    delete_2d_array(array_2d_nodes);
    delete_2d_array(array_2d_random);
    delete_2d_array(array_2d_func_in_rand);
    delete_2d_array(array_2d_func_in_nodes);
    delete_2d_array(global_matrix, Task_const::M);
    delete[] global_b_vector;
    delete_2d_array(coefficients);
    delete_2d_array(arr_approx_in_rand_points);
    delete_2d_array(array_2d_h_100_points);
    delete_2d_array(arr_approx_in_h_100_points);
    delete_2d_array(array_func_in_h_100_points);
    delete_2d_array(array_M_viz_points);
    delete_2d_array(arr_approx_in_M_viz_points);
    delete_2d_array(arr_2d_func_in_M_viz_points);

    return 0;
}