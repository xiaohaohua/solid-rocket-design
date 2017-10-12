# Solid-rocket-design
course design of the solild rocket engine
## How to use it
* install python and its libary such as numpy,scipy,matplotlib
if you are a new programmer,I recommond you to install [Anacond3](https://www.anaconda.com/download/),remember to use python3,because the code is buliding on python3
* open spider(or you can use other commander if you are familiar with python)
* change the parametres in originals.data,please pay attention to the unit of the variables.after several trials,you can find the results
at process.data,and the change of buring areas by e will be saved as graph1
## about the code
* The first part of the code is heat calculation.I suggest you to control e1 near by 0.25*D,therefore you can easily make sure y1 at the range of 0.8 and 1.2.My suggestion is to control y1 at 0.8,because the loading cofficient is limited to 0.75 to 0.85.If you get the wrong information,for I set several triggers,once your parametres is not fitable,it will be triggered.The best way to solve this problem is to change your parametres.
* The secod part and the third part is to select the design.There are several vertifies.The first is progressive cofficient,it is limited 
to less than 1.2.The second one is ![](https://github.com/graceyangfan/solid-rocket-design/raw/master/vertify2.png),for we have used arcsin 
to calculate the circumference length(s/l).The third is ![](https://github.com/graceyangfan/solid-rocket-design/raw/master/vertify3.png).
The forth is loading cofficient is limited to 0.75 to 0.85.The fifth is the left grain is limited to less than 0.05.The fininal one is the 
length of the grain.It is limited to 1.3 to 1.8.
* The finial part is to plot the graph of burning area change by e.
## plot trajectory curves
* first,you should open the original2.data to change the parametres,input D,Iz,F,pc,k,M,Tfc,density,rspeedfo(the burning speed at 6.9Mpa),n(pressure exponent),r1,tdata,L,erod.you can find their physical meaning in original2.data.
* At the first time you can set erod=1,which means you consider about the erod procession,the data is saved as pe.data.
* Then you can set erod=0,which means you do not consider about the erod procession,the data is saved as pnoe.data.
* Finally,start plot trajectory curves.py,you can get the curves.

