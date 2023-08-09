import customtkinter as ctk
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from glob import glob

#Settings
data_plotting_folder = "scripts/data_for_plotting"


# Getting info for menus
methods = [t.split("/")[-1] for t in glob(data_plotting_folder + "/*")]

exps = glob(data_plotting_folder + "/" + methods[0] + "/*")
experiments = set()
window_sizes = set()
for t in exps:
    if t.split("_")[-1].split(".")[0] == "gt":
        continue
    window_sizes.add(t.split("_")[-1].split(".")[0])
    experiments.add("_".join(t.split("/")[-1].split("_")[:-1]))
experiments = list(experiments)
experiments.sort()
window_sizes = [int(i) for i in list(window_sizes)]
window_sizes.sort()
window_sizes = [str(i) for i in list(window_sizes)]

# generate root
root = ctk.CTk()
root.geometry("1080x640")


def redraw(v):
    plt.clf()
    filename = data_plotting_folder + "/" + methodmenu.get() + "/" + filemenu.get() + "_" + windowsizemenu.get() + ".npy"
    filename_gt = data_plotting_folder + "/" + methodmenu.get() + "/" + filemenu.get() + "_" + windowsizemenu.get() + "_gt.npy"
    dets = np.load(filename)
    gt = np.load(filename_gt)
    res = (dets - gt)[np.tril_indices(dets.shape[0], k=-1)].flatten()

    cutoff = thresholdslider.get()

    inliers = res[np.abs(res) < cutoff]

    plt.hist(inliers,40)
    plt.xlabel("residual (samples)")
    plt.title(f'inlier_ratio = {inliers.size/res.size:.2f}')
    threshlabel.configure(text = f'{thresholdslider.get():.0f}')
    canvas.draw()




# generate the figure and plot object which will be linked to the root element
fig, ax = plt.subplots()
fig.set_size_inches(5,5)
#ax.hist(x,30)
#ax.axis("off")
#fig.subplots_adjust(left=0, right=1, bottom=0, top=1, wspace=0, hspace=0)

canvas = FigureCanvasTkAgg(fig,master=root)
canvas.draw()
canvas.get_tk_widget().place(anchor=ctk.E,relx=0.98, rely=0.5)

filemenu = ctk.CTkOptionMenu(master = root, values=experiments,
                                         command=redraw)
filemenu.place(relx=0.1,rely= 0.2)

methodmenu = ctk.CTkOptionMenu(master = root, values=methods,
                                         command=redraw)
methodmenu.place(relx=0.1,rely= 0.3)

windowsizemenu = ctk.CTkOptionMenu(master = root, values=window_sizes,
                                         command=redraw)
windowsizemenu.place(relx=0.1,rely= 0.4)

thresholdslider = ctk.CTkSlider(master = root,from_=10, to=300, command=redraw, width=140)
thresholdslider.place(relx=0.1,rely= 0.5)
threshlabel  = ctk.CTkLabel(master = root, text = f'{thresholdslider.get():.0f}')
threshlabel .place(relx=0.1,rely = 0.55)
redraw(1)








# Button = ctk.CTkButton(master = root,
#                        width=50,
#                        height=50,
#                        text="Gen",
#                        command=redraw)
# Button.place(relx=0.05,rely=0.425)

# slider = ctk.CTkSlider(master = root,from_=0, to=100, command=redraw)




#slider.place(relx=0.01,rely =0.1)
# initiate the window
root.mainloop()