import tkinter as tk
from tkinter import *
from tkinter import filedialog
import re
import os
column_name = ''.join("S.No" + '\t' + "Key" + '\t' +  "Start" + '\t' + "End" + '\t'  + "Identifier"  + '\t' + "Sequence")
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
        self.fasta = test2.read()

    def reset(self):
        root.destroy()

    def save(self):
        text_file = open("output.txt", "w")
        text_file.write(text_box.get(1.0, tk.END))
        text_file.close()

    def intersection(self):
        signal_adjusted_fasta = []
        count = 0
        for signal in self.file_one:

            for x in self.fasta.strip().split('>'):
                content = x.split('\n')
                # Since after splitting by '>' first element will be blank so len(element) >1 is required
                if len(content) > 1:
                    header = content[0]
                    sequence = str(''.join(content[1:]))
                    sequence = sequence.upper()
                    # print(sequence)
                    # adjusted_fasta = ''.join(header + '\t' + sequence)
                    # data.append(adjusted_fasta)
                    matched = re.findall(signal, sequence)
                    for match in re.finditer(signal, sequence):
                        count += 1
                        if len(matched) > 0:
                            # data = ''.join(str(count) + '\t' + str(','.join(matched)) + '\t' +  str(match.start()) + '\t' + str(match.end()) + '\t'  + header  + '\t' + sequence)
                            data = ''.join(str(count) + '\t' + str(','.join(matched)) + '\t' +  str(match.start()) + '\t' + str(match.end()) + '\t'  + header  + '\t' + sequence)
                            signal_adjusted_fasta.append(data)

        text_box.insert(tk.END, column_name +'\n')
        for k in signal_adjusted_fasta:
            text_box.insert(tk.END, k +'\n')
        count = 0
        for m in signal_adjusted_fasta:
            count += 1
        count_box.insert(tk.END, str(count) + " key signal peptide were found to be present in given fasta")


root = tk.Tk()
root.title(' ----- Signal Peptide  Finder Tool -----')
root.geometry("1500x800")
root.configure(background='gray89')
a = ex()
# This will create a LabelFrame
label_frame = LabelFrame(root, text = 'Signal Peptide Finder Tool', width= 300, height= 200, labelanchor= "n",font= ('Helvetica 14 bold', 40),bd= 10, background="white", foreground= "Black")
label_frame.pack(ipadx=0, ipady=0,  expand = True, fill = "x" )
#Create a Label inside LabelFrame
Label(label_frame, text= "Instructions!!\nUpload peptide signal file (Format should be regular expression pattern) and fasta of your interest.\n Then click Explore and Click Save. \n Must press Reset at the end, which will close the app and need to be started again", font=('Helvetica 15  underline', 15), foreground= "black").pack(pady= 20)

button1 = tk.Button(root, text='Enter peptide signal File', command=a.uploadaction1, relief=GROOVE)
button2 = tk.Button(root, text='Enter Fasta File', command=a.uploadaction2, relief=GROOVE)
inter = tk.Button(root, text="Explore", command=a.intersection, relief=RAISED)
save = tk.Button(root, text="Save", command=a.save, relief=RAISED)
reset = tk.Button(root, text="Reset", command=a.reset, relief=RAISED)

button1.pack()
pathlabel1 = Label(root)
pathlabel1.pack()
button2.pack()
pathlabel2 = Label(root)
pathlabel2.pack()
inter.pack()
save.pack()
reset.pack()

count_box = tk.Text(height =1, width = 200, bg = "light yellow", font=('Arial', 13))
count_box.pack()
text_box = tk.Text(height =30, width = 200, bg = "lavenderblush1", font=('Arial', 12))
text_box.pack()

root.mainloop()
