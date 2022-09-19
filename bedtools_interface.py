import tkinter as tk
from tkinter import *
from tkinter import filedialog
class ex:
    def uploadaction1(self):
        filename1 = filedialog.askopenfilename()
        file1label = filename1.split("/")
        file1label = file1label[-1]
        pathlabel1.config(text=file1label, fg="blue")
        test1 = open(filename1)
        self.file_one = test1.read().strip().split("\n")

    def uploadaction2(self):
        filename2 = filedialog.askopenfilename()
        pathlabel2.config(text=filename2, fg="navy")
        file2label = filename2.split("/")
        file2label = file2label[-1]
        pathlabel2.config(text=file2label, fg="darkred")
        test2 = open(filename2)
        self.file_two = test2.read().strip().split("\n")

    def intersection(self):
        output_data = []

        #         print(self.file_one)

        for i in self.file_one:
            i = i.split("\t")
            #             print("splitted file",i)

            for j in self.file_two:
                j = j.split("\t")

                if i[0] == j[0]:

                    if int(i[1]) <= int(j[1]) and int(i[2]) >= int(j[1]):
                        #                         print("first:","\t".join(i + j))
                        output_data.append("\t\t".join(i + j))
                    elif int(i[1]) >= int(j[1]) and int(i[2]) <= int(j[2]):
                        output_data.append("\t\t".join(i + j))
                    #                         print("second:","\t".join(i + j))
                    elif int(i[1]) >= int(j[1]) and int(i[1]) <= int(j[2]):
                        output_data.append("\t\t".join(i + j))
                    #                         print("3:","\t".join(i + j))
                    elif int(i[1]) <= int(j[1]) and int(i[2]) >= int(j[2]):
                        output_data.append("\t\t".join(i + j))
                    #                         print("4","\t".join(i + j))
                    elif int(j[1]) >= int(i[1]) and int(j[1]) <= int(i[2]):
                        output_data.append("\t\t".join(i + j))
                    #                         print("5","\t".join(i + j))
                    elif int(j[1]) <= int(i[1]) and int(j[2]) >= int(i[1]):
                        output_data.append("\t\t".join(i + j))
        #                         print("6","\t".join(i + j))
        count = 0
        for k in output_data:
            text_box.insert(tk.END, k +'\n')
        count = 0
        for m in output_data:
            count += 1
        count_box.insert(tk.END, "Total elements overlapped between two files are:   "+ str(count))



root = tk.Tk()
root.title('BED Intersector')
root.geometry("1500x700")
root.configure(background='gray89')
a = ex()
# This will create a LabelFrame
label_frame = LabelFrame(root, text = 'BED File Intersection Tool', width= 200, height= 200, labelanchor= "n",font= ('Helvetica 14 bold', 40),bd= 10, background="white", foreground= "Black")
label_frame.pack(ipadx=0, ipady=0,  expand = True, fill = "x" )
#Create a Label inside LabelFrame
Label(label_frame, text= "Instructions!! Upload two bed files of your interest and click intersect", font=('Helvetica 15 bold', 20), foreground= "black").pack(pady= 20)

button1 = tk.Button(root, text='Enter Bed File 1', command=a.uploadaction1, relief=GROOVE)
button2 = tk.Button(root, text='Enter Bed File 2', command=a.uploadaction2, relief=GROOVE)
inter = tk.Button(root, text="Intersect", command=a.intersection, relief=RAISED)
button1.pack()
pathlabel1 = Label(root)
pathlabel1.pack()
button2.pack()
pathlabel2 = Label(root)
pathlabel2.pack()
inter.pack()
count_box = tk.Text(height =1, width = 160, bg = "light yellow", font=('Arial', 15))
count_box.pack()
text_box = tk.Text(height =30, width = 180, bg = "lavenderblush1", font=('Arial', 12))
text_box.pack()

root.mainloop()
