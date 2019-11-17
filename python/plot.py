import matplotlib.pyplot as plt

def get_data(file_name):
    content = []
    with open(file_name) as f:
        content = f.readlines() # read file line by line
    content = [x.strip() for x in content]
    return content

def get_info(content, i, j, k):
    sizes = []
    algo = []
    ward = []
    for e in range(len(content)):
        c = content[e].split()
        if e % 3 == 0:
            sizes.append(int(c[i]))
        if e % 3 == 1:
            algo.append(float(c[j]))
        if e % 3 == 2:
            ward.append(float(c[k]))
    return sizes, algo, ward

def get_info_perf(content, i, j):
    x = []
    y = []
    for e in range(len(content)):
        c = content[e].split()
        x.append(int(c[i]))
        y.append(float(c[j]))
    return x, y

# content = get_data('result100_25_128.txt')
# sizes, algo, ward = get_info(content)
#
# plt.plot(sizes, algo, label = 'Ward-approx')
# plt.plot(sizes, ward, label = 'Ward')
# plt.title('Accuracy of Ward-approx (e = 10, T = 25, L = 128) vs Ward')
# plt.legend(loc = 'best')
# plt.xlabel('Number of Points')
# plt.ylabel('Accuracy NMI')
# plt.show()

# content = get_data('result_epsilons.txt')
# eps, alg_acc, ward_acc = get_info(content, 1, 1, 1)
# plt.plot(eps, alg_acc, label = 'Ward-approx')
# plt.title('The effect of epsilon on the Accuracy')
# plt.legend(loc = 'best')
# plt.xlabel('Epsilons')
# plt.ylabel('Accuracy NMI')
# plt.show()


def param_vs_perf(perf_file, i, j, title, x_label, y_label):
    content = get_data(perf_file)
    x, y = get_info_perf(content, i, j)
    plt.plot(x, y, label = 'Ward-approx')
    plt.title(title)
    plt.legend(loc = 'best')
    plt.xlabel(x_label)
    plt.ylabel(y_label)
    plt.show()



#param_vs_perf("leaves_newsgroup_perfs.txt", 2, 5, 'The effect of the number of leaves on the performance', 'Number of Leaves', 'Performance in sec')
def params_vs_acc(acc_file, i, j, k, title, x_label, y_label):
    content = get_data(acc_file)
    x, alg_acc, ward_acc = get_info(content, i, j, k)
    plt.plot(x, alg_acc, label = 'Ward-approx')
    plt.title(title)
    plt.legend(loc = 'best')
    plt.xlabel(x_label)
    plt.ylabel(y_label)
    plt.show()

params_vs_acc('result_leaves.txt', 3, 1, 1, 'The effect of the number of leaves on the accuracy', 'the Number of Leaves', 'the Accuracy')
