#!/usr/bin/python
# TRAX (c) 2012, Lomonosov Moscow State University
# http://erg.biophys.msu.ru/wordpress/archives/330

import matplotlib
matplotlib.use('TkAgg')

from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2TkAgg
from matplotlib.figure import Figure
from matplotlib.colors import rgb2hex,colorConverter

from Tkinter import *
import tkMessageBox, tkSimpleDialog, tkFont
from tkFileDialog import askopenfilename, asksaveasfilename

import numpy as np
from scipy.integrate import odeint
from copy import copy

import re
import sys
import traceback
import math

class SODE:
  def __init__(self, text):
    self.var=[]
    self.var1=[]
    self.dvar=[]
    self.rhs={}
    self.drhs={}
    self.t=np.array([])
    self.kin={}
    self.kin1={}

    for s in text.split('\n'):
      _s = s.strip()
      if not len(_s):
        continue

      # skip comments
      if _s[0] == '#':
        continue

      if '=' in s:
        lhs,rhs=map(lambda s: s.replace(' ',''), s.split('='))
        if lhs[-1]=="'":
          v=lhs[:-1]
          self.drhs[v]=rhs
          if not v in self.dvar: self.dvar.append(v)
        else:
          v=lhs
          self.rhs[v]=rhs
        if not v in self.var: self.var.append(v)

    for v in self.var:
      rhv = re.sub("[-+/*()]"," ",self.rhs[v]).split()
      for v1 in rhv:
        if v1 in self.dvar+['t']:
          self.var1.append(v)
          break

    self.par=[]
    for v in self.var:
      if not v in (self.dvar+self.var1): self.par.append(v)


    self.parval={}
    for p in self.par:
      try:
        self.parval[p] = eval(self.rhs[p], dict(math.__dict__.items() + self.parval.items()))
      except StandardError as err:
        raise sodeError, "Error parsing following expression:\n`%s = %s`\n%s"%(p, self.rhs[p], err.message)

    self.f={}
    # add various math functions to the namespace of the expression for the equation RHS
    nmsps = dict(math.__dict__.items() + self.parval.items())
    
    for v in self.dvar:
      s=self.drhs[v]
      i=0
      for _v in self.dvar:
        nsub=1
        while nsub: s,nsub=re.subn("(^|[-+*/(])%s([-+*/)]|$)" % _v, "\\1x[%d]\\2" % i, s)
        i+=1
#      print v, s
      if not v in self.rhs.keys(): raise sodeError, "Initial value of `%s` is undefined"%v
      try:
        self.f[v]=eval("lambda x,t: %s"%s, nmsps)
      except StandardError as err:
        raise sodeError, "Error parsing following expression:\n`%s' = %s`\n%s"%(v, self.drhs[v], err.message)
    for v in self.dvar:
      try:
        # make a test call to a function to workaround odeint() behavior, which doesn't raise exceptions on errors.
        tmp=self.f[v](map(lambda v: eval(self.rhs[v], dict(math.__dict__.items() + self.parval.items())), self.dvar),0)
      except Exception as err:
        raise sodeError, "Error evaluating the following expression:\n `%s' = %s`\n%s"%(v, self.drhs[v], str(err))

  def getdx(self):
    if self.parval['t']<0: sig=-1
    else: sig = 1
    return lambda _x,_t: sig*np.array(map(lambda v,x=_x,t=_t: self.f[v](x,t), self.dvar))

  def integrate(self):
    atol=1.49012e-8
    rtol=1.49012e-8
    if len(self.dvar)==0: raise sodeError, "No equations defined."
    if 't' not in self.parval.keys(): raise sodeError, "Integration time limit 't' is undefined."
    #TODO: implement adaptive grid
    t = np.linspace(0,self.parval['t'], 10) # initial 
    t.sort() # sort in case t value was negative
    dx = self.getdx()
    parvals = map(lambda v: eval(self.rhs[v], dict(math.__dict__.items() + self.parval.items())), self.dvar)
    x = odeint(dx, parvals, t,atol=atol,rtol=rtol) # integrate ODEs
    if np.isnan(x).sum(): raise NameError, "Integration failed" # check for NANs in output
    D0 = 0 # curve length from previous iteration
    refine = True # flag to refine grid
    npadd = 2
    print "\n\nStart grid refinement"
    iref=0 # iteration counter
    while t.shape[0]<10000 and refine: # limit number of grid points
      iref += 1
      ds = map(lambda i: np.linalg.norm(x[i+1,:] - x[i,:]), range(x.shape[0]-1)) # curve segments lengths
      D = sum(ds) # curve length
      if (D == 0): break # stationary point?
      dmax = max(ds) # max segment length
#      C = x.mean(0) # curve center
#      rs = map(lambda i: np.linalg.norm(x[i,:] - C), range(x.shape[0])) # curve points distance from center
#      S = max(rs) # curve radius
      print "\nGrid points count: %d"%t.shape[0]
      print "Track length: %f"%D
      print "Track length change: %.1f%%"%((1-D0/D)*100)
      print "Max segment length: %f"%max(ds)
#      print "Track size: %f"%S
      refine = False
      for i in xrange(0,t.shape[0]-1):
        a=x[i+1,:]-x[i,:]
        b=dx(x[i,:], t[i])*(t[i+1]-t[i])
        d = np.linalg.norm(a-b)/np.linalg.norm(a)
        if d > 1e-1 and np.abs((a / (atol + rtol*x[i]))).max() > 1:
          t.resize((t.shape[0]+npadd))
          t[-npadd:] = np.linspace(t[i],t[i+1],npadd+2)[1:-1]
#          t.resize((t.shape[0]+1))
#          t[-1] = (t[i]+t[i+1])/2
#          print "refine at %f (%e,%f)"%(t[-1], np.linalg.norm(b-a), d)
          refine = True
      if not refine and abs((D-D0)/D) > 5e-2:
        l=t.shape[0]
        t.resize((2*l-1))
        t[-l+1:] = (t[1:l] +t[:l-1])/2
        refine = True
      D0 = D
      t.sort()
      x = odeint(dx, parvals, t, atol=atol,rtol=rtol)
      print "New grid points count: %d"%t.shape[0]
    if not refine: print "Grid is fine",
    else: print "Grid points limit reached",
    print " after %d iterations"%iref
    self.t = t
    self.tracks = x
    self.tracks2dict()

  def tracks2dict(self):
    i=0
    for v in self.dvar: 
      self.kin[v] = self.tracks[:,i]
      i+=1
    for v in self.var1:
      self.kin[v] = eval(self.rhs[v],dict(self.kin.items()+self.parval.items() + [('t',self.t)]))
    self.dvar = self.dvar+self.var1 # TODO: NOT SAFE! second call to integrate() will fail.
      

  def discrete(self):
    if len(self.dvar)==0: raise sodeError, "No equations defined."
    if 't' not in self.parval.keys(): raise sodeError, "Integration time limit 't' is undefined."
    self.t=np.arange(self.parval['t'])
    self.tracks = np.zeros((self.t.shape[0], len(self.dvar)))
    self.tracks[0,:] = map(lambda v: eval(self.rhs[v], self.parval), self.dvar)
    dx = self.getdx()
    for i in self.t[1:]:
      self.tracks[i,:] = dx(self.tracks[i-1,:], i)
    self.tracks2dict()

  def __getitem__(self, n):
    if n=='t': return copy(self.t)
    else:
      if n in self.dvar: return self.kin[n]
      raise KeyError

class sodeError(StandardError):
  pass

class Trax_GUI:
  def __init__(self, master):
    self.discrete = False
    self.mpl_ver_maj=int(matplotlib.__version__.split('.')[0])
    self.master = master
    self.font_size=12
    self.min_height = int(390*self.font_size/10.)
    self.line_height = int(27*self.font_size/10.)
    self.height = 0.95*self.master.winfo_screenheight()
    if self.height > 2*self.min_height: self.height=2*self.min_height
    if self.height <= 700: self.font_size=10
    master.wm_title("TraX")
    fnt=tkFont.Font(font=("Helvetica",self.font_size,NORMAL))
    master.option_add("*Font", fnt)

    # Buttons
    fr_but = Frame(master)
    fr_but.pack(side=TOP, anchor='nw')

    Button(fr_but, text="Open", command=self.load).pack(side=LEFT)
    Button(fr_but, text="Save", command=self.save).pack(side=LEFT)

    fr_fig_sel = Frame(master)
    fr_fig_sel.pack(side=TOP, anchor='center')
    self.fr_fig_sel=fr_fig_sel
    self.figlist=[]
    Button(fr_fig_sel, command = self.newfig, text="New fig", bg="cyan").pack(side=LEFT)

    fr_fig_mnu = Frame(master)
    fr_fig_mnu.pack(side=TOP, anchor='nw')


    # Figure
    fr_fig = Frame(master,padx=5, pady=5)
    fr_fig.pack(side=TOP)

    fr_fig1=Frame(fr_fig)
    fr_fig1.pack(side=LEFT)
    self.f = Figure(figsize=(5.5,(self.height-self.min_height+2*self.line_height)/100.), dpi=100)
    canvas = FigureCanvasTkAgg(self.f, master=fr_fig1)
    canvas.show()
    # Figure toolbar
    toolbar = NavigationToolbar2TkAgg(canvas, fr_fig1)
    toolbar.update()
    canvas.get_tk_widget().pack(side=LEFT)
    self.ax = self.f.add_subplot(111)
    # Figure controls
    fr_plt_ctrl=Frame(fr_fig, width=150)
    fr_plt_ctrl.pack(side=LEFT, anchor='nw')
    self.fr_plt_ctrl=fr_plt_ctrl

    fr_eq = Frame(master, padx=5, pady=5)
    fr_eq.pack(side=TOP, anchor='nw')

    # Equations input 
    self.txt=Text(fr_eq, width=80, height=7)
    self.txt.pack(side=LEFT)
    self.txt.bind("<<Modified>>", self.alert_eq_mod)
    self.rst_mod_flg=False
    scrbar=Scrollbar(fr_eq, command=self.txt.yview)
    scrbar.pack(side=LEFT, fill=Y)
    self.txt["yscrollcommand"]=scrbar.set

    # RUN RABBIT RUN
    fr_run = Frame(master)
    fr_run.pack(side=TOP)
    self.but_run = Button(fr_run, text="RUN", bg="green", fg="black", width=5, command=self.cauchy)
    self.but_run.pack(side=LEFT)
    self.discrete = IntVar()
    Checkbutton(fr_run,text="Discrete",variable=self.discrete).pack(side=LEFT)


    self.ax.figure.canvas.mpl_connect("button_press_event", self.get_init_vals)
    self.do_get_init_vals = False
    self.but_get_init_vals=Button(fr_fig_mnu, text="Point and Run", command=self.start_get_init_vals)
    self.but_get_init_vals.pack(side=LEFT)

    Label(fr_fig_mnu, text="Xmin").pack(side=LEFT)
    self.en_xmin = Entry(fr_fig_mnu, width=4)
    self.en_xmin.pack(side=LEFT)
    self.en_xmin.lim_id = 0

    Label(fr_fig_mnu, text="Xmax").pack(side=LEFT)
    self.en_xmax = Entry(fr_fig_mnu, width=4)
    self.en_xmax.pack(side=LEFT)
    self.en_xmax.lim_id = 1

    Label(fr_fig_mnu, text="Ymin").pack(side=LEFT)
    self.en_ymin = Entry(fr_fig_mnu, width=4)
    self.en_ymin.pack(side=LEFT)
    self.en_ymin.lim_id = 2

    Label(fr_fig_mnu, text="Ymax").pack(side=LEFT)
    self.en_ymax = Entry(fr_fig_mnu, width=4)
    self.en_ymax.pack(side=LEFT)
    self.en_ymax.lim_id = 3

    self.en_lim=(self.en_xmin,self.en_xmax,self.en_ymin,self.en_ymax)

    self.sode = SODE("")
    self.newfig()


  def newfig(self):
    self.figlist.append(Fig(Button(self.fr_fig_sel),self))
    self.figlist[-1].activate()

  def alert_eq_mod(self, event):
    # alert SODE modification by user
    if self.rst_mod_flg: return
    self.rst_mod_flg=True
    self.txt.tk.call(self.txt._w, "edit", "modified", 0)
    self.rst_mod_flg=False
    self.but_run['bg']="yellow"

  def quit(self):
    if tkMessageBox.askyesno("Quit", "Really?"): self.master.quit()

  def cauchy(self):
    # parse and solve newly defined SODE
    try:
      self.sode=SODE(self.txt.get(1.0,END))
      if self.discrete.get(): self.sode.discrete()
      else: self.sode.integrate()
      
    except sodeError as err:
      tkMessageBox.showerror("ERROR", err.message)
    except StandardError as err:
      tkMessageBox.showerror("ERROR", "An error occured while parsing your input:\n%s\n--\nException arguments:\n %s\n%s"%(err.message, str(err.args), traceback.format_exc()))
    else:
      # update figures data
      for f in self.figlist:
        f.sode=self.sode
        # update list of variables in current control frame
        if f.active: f.activate() 
      self.but_run["bg"]="green" # make "RUN" button green

  def start_get_init_vals(self):
    if not self.do_get_init_vals:
      self.but_get_init_vals['bg']="yellow"
      self.do_get_init_vals = True
    else:
      self.do_get_init_vals=False
      self.but_get_init_vals['bg']="grey"


  def get_init_vals(self, event):
    if self.do_get_init_vals:
      for f in self.figlist:
        if f.active: break
      if len(f.curves):
        for n,v in ((f.curves[-1].x,event.xdata),(f.curves[-1].y,event.ydata)):
          if not n == 't':
            self.txt.insert(END, "\n%s = %f\n"%(n,v))
        self.cauchy()
        f.add_curve(f.curves[-1].x,f.curves[-1].y)


  def load(self):
    fn=askopenfilename()
    if fn:
      f=open(fn)
      self.txt.delete(1.0, END)
      self.txt.insert(END, f.read())
      f.close()

  def save(self):
    fn=asksaveasfilename()
    if fn:
      f=open(fn, 'w')
      print >>f, self.txt.get(1.0, END)
      f.close()


class Fig:
  def __init__(self, button, trax):
    self.trax=trax

    self.fl=trax.figlist
    self.sode=trax.sode
    self.ax=trax.ax
    self.fr=trax.fr_plt_ctrl
    self.but=button
    self.but["text"]="Fig"
    self.but["command"]=self.activate
    self.but.pack(side=LEFT)
    self.curves=[]

    self.lim=[0,0,0,0]
    self.apply_lim=[0,0,0,0]
    self.en_lim=trax.en_lim # list of tk.Entries to input axes limits
    self.x,self.y='t','t'
    self.activate()

  def test_lim(self, event):
    if self.active:
      if len(event.widget.get().strip()):
        self.lim[event.widget.lim_id]=float(event.widget.get())
        self.apply_lim[event.widget.lim_id]=1
#        print "apply lim %d"%event.widget.lim_id
        self.scale()
      else:
        self.apply_lim[event.widget.lim_id]=0
#        print "disable lim %d"%event.widget.lim_id
  

  def add_curve(self,x,y):
    self.x,self.y=x,y
    self.curves.append(Curve(self.sode, x, y, self))
    self.scale()

  def activate(self):
    # bind events from limits entries to own methods
    for i in xrange(4):
      self.en_lim[i].bind("<FocusOut>", self.test_lim)
    # highlight own window button
    self.active=True
    for f in self.fl:
      if f is self:
        self.but["bg"]="yellow"
      else:
        f.but["bg"]="grey"
        f.active=False
    # clear and repopulate plot control frame
    for item in self.fr.pack_slaves():
      item.pack_forget()
      item.destroy()
    Button(self.fr, text="Delete figure",bg="red",fg="white",command=self.destroy).pack(side=TOP)
    f=Frame(self.fr)#, bg="magenta", padx=3, pady=3)
    f.pack(side=TOP)
    if not self.x in self.sode.dvar: self.x='t' 
    if not self.y in self.sode.dvar: self.y='t' 
    v1=Variable()
    v1.set(self.y)
    v2=Variable()
    v2.set(self.x)
    OptionMenu(f, v1, 't',*self.sode.dvar).pack(side=LEFT)
    Label(f, text=" (").pack(side=LEFT)
    OptionMenu(f, v2, 't',*self.sode.dvar).pack(side=LEFT)
    Label(f, text=") ").pack(side=LEFT)
    Button(f, text="Add", command=lambda: self.add_curve(v2.get(),v1.get())).pack(side=LEFT)

    f=Frame(self.fr)#, bg="magenta", padx=3, pady=3)
    f.pack(side=TOP)
    self.cnvs=Canvas(f, bg="grey70", width=120, height=self.trax.height-self.trax.min_height)
    self.cnvs.pack(side=LEFT)
    scrbr=Scrollbar(f, command=self.cnvs.yview)
    scrbr.pack(side=LEFT, fill=Y)
    self.cnvs["yscrollcommand"]=scrbr.set
    self.fr_curves = Frame(self.cnvs, bg="grey70")
    self.cnvs.create_window((5,5), window=self.fr_curves, anchor='nw')
    Button(self.fr, text="Clear figure", bg="white", command=self.clear).pack(side=TOP)
#    self.fr_curves.pack()

    # clear plot
    self.ax.cla()
    self.ax.grid()
    # restore curves and their controls
    for c in self.curves:
      c.show()
    # restore axes limits
    for i in xrange(4):
      self.en_lim[i].delete(0,END)
      if self.apply_lim[i]:
        self.en_lim[i].insert(0, "%f"%self.lim[i])
    # scale and redraw figure
    self.scale()

  def destroy(self):
    # if there are other figures
    if len(self.fl)>1:
    # remove own window button
      for i in xrange(len(self.fl)):
        if self.fl[i] is self:
          del self.fl[i]
          if i: self.fl[i-1].activate()
          else: self.fl[0].activate()
          break
      self.but.pack_forget()
      self.but.destroy()

  def scale(self):
    # scale figure apply limits and redraw
    if self.trax.mpl_ver_maj: self.ax.autoscale()
    lim=[0,0,0,0]
    lim[:2] = self.ax.get_xlim()
    lim[2:] = self.ax.get_ylim()
    for i in xrange(4):
      if self.apply_lim[i]:
        lim[i]=self.lim[i]
    self.ax.set_xlim(lim[:2])
    self.ax.set_ylim(lim[2:])
    self.ax.figure.canvas.draw()

  def clear(self):
    for c in copy(self.curves):
      c.destroy()


class Curve:
  def __init__(self, sode, x, y, fig):
    self.x,self.y=x,y
    self.sode=sode
    self.fig = fig
    self.show()

  def show(self):
    self.fr1=Frame(self.fig.fr_curves, padx=3, pady=3)
    self.fr1.pack(side=TOP)
    l=Label(self.fr1, text="%s(%s)"%(self.y, self.x)); l.pack(side=LEFT)
    Button(self.fr1, text="Del", command=self.destroy).pack(side=LEFT)
#    self.line, = self.fig.ax.plot(self.sode[self.x], self.sode[self.y],'.-',label="%s(%s)"%(self.y, self.x))
    self.line, = self.fig.ax.plot(self.sode[self.x], self.sode[self.y],label="%s(%s)"%(self.y, self.x))
    self.fr1["bg"] = rgb2hex(colorConverter.to_rgb(self.line.get_color()))
#    self.ax.legend()
#    self.ax.figure.canvas.draw()
    self.fig.cnvs.update_idletasks()
    self.fig.cnvs["scrollregion"]=self.fig.cnvs.bbox(ALL)

  def destroy(self):
    self.fig.ax.lines.remove(self.line)
#    self.ax.legend()
    self.fig.ax.figure.canvas.draw()
    i=0
    for i in xrange(len(self.fig.curves)):
      if self.fig.curves[i] is self:
        del self.fig.curves[i]
        break
    self.fr1.pack_forget()
    self.fr1.destroy()
    self.fig.cnvs["scrollregion"]=self.fig.cnvs.bbox(ALL)

if __name__ == "__main__":
  root = Tk()
  app = Trax_GUI(root)
  root.mainloop()
