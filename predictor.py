import tkinter as tk
from tkinter import *
from tkinter import filedialog
import re
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


        for k in signal_adjusted_fasta:
            text_box.insert(tk.END, k +'\n\n\n')
        count = 0
        for m in signal_adjusted_fasta:
            count += 1
        count_box.insert(tk.END, "Total elements overlapped between two files are:   "+ str(count))


root = tk.Tk()
root.title('Peptide Signal Finder')
root.geometry("1500x700")
root.configure(background='gray89')
a = ex()
# This will create a LabelFrame
label_frame = LabelFrame(root, text = 'Peptide Signal Predictor Tool', width= 300, height= 200, labelanchor= "n",font= ('Helvetica 14 bold', 40),bd= 10, background="white", foreground= "Black")
label_frame.pack(ipadx=0, ipady=0,  expand = True, fill = "x" )
#Create a Label inside LabelFrame
Label(label_frame, text= "Instructions!! Upload peptide key (Format: follow regular expression pattern) and fasta of your interest and click Explore", font=('Helvetica 15 bold', 20), foreground= "black").pack(pady= 20)

button1 = tk.Button(root, text='Enter peptide signal File', command=a.uploadaction1, relief=GROOVE)
button2 = tk.Button(root, text='Enter Fasta File', command=a.uploadaction2, relief=GROOVE)
inter = tk.Button(root, text="Explore", command=a.intersection, relief=RAISED)
button1.pack()
pathlabel1 = Label(root)
pathlabel1.pack()
button2.pack()
pathlabel2 = Label(root)
pathlabel2.pack()
inter.pack()
count_box = tk.Text(height =1, width = 200, bg = "light yellow", font=('Arial', 15))
count_box.pack()
text_box = tk.Text(height =30, width = 200, bg = "lavenderblush1", font=('Arial', 12))
text_box.pack()

root.mainloop()