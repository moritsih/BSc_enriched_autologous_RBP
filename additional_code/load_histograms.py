


import numpy as np


def my_hist(data, borders = "standard", bin_number = "standard", window_average = 5, gauss = True):
    if window_average%2 !=1:
        raise ValueError('window-size must be uneven!')

    data = [float(x) for x in data if not np.isnan(x)]
    sigma = np.std(data)
    mean = np.mean(data)
    median = np.median(data)
    n = len(data)

    if borders == "standard":
        borders = [mean-5*sigma, mean+5*sigma]

    if borders == "full":
        borders = [min(data), max(data)]

    if bin_number == "standard":
        bin_size = 3.49*sigma*n**(-1/3)
        bin_number = int((borders[1]-borders[0])/bin_size) +1

    x = np.linspace(borders[0],borders[-1],bin_number+2)[1:-1]
    x_smooth = [x[0]]

    for i in range(len(x)-1):
        x_smooth.append((x[i]+x[i+1])/2)
        x_smooth.append(x[i+1])

    y, y_smooth = [0]*bin_number, [0]*bin_number

    #fill
    for number in data:
        if number<borders[0]:
            y[0]+=1
        elif number>=borders[1]:
            y[bin_number-1]+=1
        else:
            y[int((number-borders[0])//((borders[1]-borders[0])/bin_number))] += 1


    #interpolation
    y_smooth = []
    for i in range(len(y)-1):
        if i==0 or i==len(y)-2:
            y_smooth.append(y[i])
            y_smooth.append(np.mean(y[i:i+2]))
        else:
            y_smooth.append(y[i])
            y_smooth.append(np.mean(y[i-1:i+3]))
    y_smooth.append(y[-1])


    #smoothing
    for i in range(len(y_smooth)):
        if i < window_average//2 or i>=len(y_smooth)-window_average//2:
            y_smooth[i] = y_smooth[i]/window_average*(window_average//2+1)
        else: 
            y_smooth[i] = np.mean(y_smooth[i-window_average//2:i+window_average//2+1])

    #area=1
    y_smooth, y = [yi/sum(y_smooth)*len(y_smooth)/len(y) for yi in y_smooth], [yi/sum(y) for yi in y]

    width=(borders[1]-borders[0])/bin_number*0.95
    if gauss:
        #gauss_fit
        e = 2.71828
        a = max(y)

        #improve a
        y_arr = np.array(y)
        diff = 10**10
        for a_improved in np.arange(a*0.3,a*1.3,a/1000):
            y2 = [a_improved*e**(-(((ex-mean)**2)/2/(sigma**2))) for ex in x]
            y_calc_arr = np.array(y2)
            if sum(abs(y_arr-y_calc_arr))<diff:
                best_a = a_improved
                diff = sum(abs(y_arr-y_calc_arr))

        gauss_x = [x for x in np.arange(borders[0],borders[1],(borders[1]-borders[0])/1000)]
        gauss_y = [best_a*e**(-(((x-mean)**2)/2/(sigma**2))) for x in gauss_x]
        return({'x':x, 'y':y,'x_smooth':x_smooth,'y_smooth':y_smooth, 'gauss_x':gauss_x, 'gauss_y':gauss_y, 'width':width, 'mean':mean, 'median':median, 'sigma':sigma, 'bin_number':bin_number, 'borders':borders})
    return({'x':x, 'y':y,'x_smooth':x_smooth,'y_smooth':y_smooth, 'width':width, 'mean':mean, 'median':median, 'sigma':sigma, 'bin_number':bin_number, 'borders':borders})



def vertical_hist(data, xpos, x_width, borders = "full", bin_number = "standard", window_average = 5):
    if window_average%2 !=1:
        raise ValueError('window-size must be uneven!')

    if not isinstance(data[0], list):
        data = [float(x) for x in data if not np.isnan(x)]
    else:
        data = [float(x) for sub in data for x in sub if not np.isnan(x)]

    sigma = np.std(data)
    mean = np.mean(data)
    median = np.median(data)
    n = len(data)

    if borders == "standard":
        borders = [mean-5*sigma, mean+5*sigma]

    if borders == "full":
        borders = [min(data), max(data)]

    if bin_number == "standard":
        bin_size = 3.49*sigma*n**(-1/3)
        bin_number = int((borders[1]-borders[0])/bin_size)+1

    #edges = np.linspace(borders,borders,bin_number+1)
    x = np.linspace(borders[0],borders[-1],bin_number+2)[1:-1]
    x_smooth = [x[0]]

    for i in range(len(x)-1):
        x_smooth.append((x[i]+x[i+1])/2)
        x_smooth.append(x[i+1])

    y, y_smooth = [0]*bin_number, [0]*bin_number

    #fill
    for number in data:
        if number<borders[0]:
            y[0]+=1
        elif number>=borders[1]:
            y[bin_number-1]+=1
        else:
            y[int((number-borders[0])//((borders[1]-borders[0])/bin_number))] += 1


    #interpolation
    y_smooth = []
    for i in range(len(y)-1):
        if i==0 or i==len(y)-2:
            y_smooth.append(y[i])
            y_smooth.append(np.mean(y[i:i+2]))
        else:
            y_smooth.append(y[i])
            y_smooth.append(np.mean(y[i-1:i+3]))
    y_smooth.append(y[-1])


    #smoothing
    for i in range(len(y_smooth)):
        if i < window_average//2 or i>=len(y_smooth)-window_average//2:
            y_smooth[i] = y_smooth[i]/window_average*(window_average//2+1)
        else: 
            y_smooth[i] = np.mean(y_smooth[i-window_average//2:i+window_average//2+1])

    #area=1
    y_smooth, y = [yi/sum(y_smooth)*len(y_smooth)/len(y) for yi in y_smooth], [yi/sum(y) for yi in y]


    #make_vertical
    vert_x, vert_y = [], x
    max_width = max(y)
    for x_new in y:
        vert_x.append(x_new/max_width*x_width+xpos)

    vert_x_smooth, vert_y_smooth = [], x_smooth
    max_width = max(y_smooth)
    for x_new in y_smooth:
        vert_x_smooth.append(x_new/max_width*x_width+xpos)

    return {'x':vert_x, 'y':vert_y, 'x_smooth':vert_x_smooth,'y_smooth':vert_y_smooth, 'mean':mean, 'median':median, 'sigma':sigma, 'bin_number':bin_number, 'borders':borders}








# data = [1.5]*2+[2.5]*1+[3.5]*4+[4.5]*8+[5.5]*15+[6.5]*11+[7]*5+[8.5]*4+[9.5]*2+[10.5]*1
# x,y,x_smooth,y_smooth,gauss_x, gauss_y, width, mean, median, sigma, bin_number = my_hist(data, borders = [1,11], bin_number = 10, window_average = 3, gauss = True)

# plt.bar(x,y, width=width)
# plt.plot(x_smooth, y_smooth, color='b')
# plt.plot(gauss_x, gauss_y, color='r')
# plt.xticks([1,2,3,4,5,6,7,8,9,10])
# plt.show()



# data = [1.5]*2+[2.5]*1+[3.5]*4+[4.5]*8+[5.5]*15+[6.5]*11+[7]*5+[8.5]*4+[9.5]*2+[10.5]*1
# vert_x, vert_y, vert_x_smooth, vert_y_smooth, mean, median = vertical_hist(data, 1, 2, borders = "full", bin_number = 10, window_average = 5)

# plt.plot(vert_x, vert_y, color='b')
# plt.plot(vert_x_smooth, vert_y_smooth, color='r')

# plt.show()
