<div align="left">
<img style="flex-grow:1; flex-shrink:1; border: 0px solid black;" src="./logo.jpg" width="210" />
</div>
<div align="left">

[LiuLab]() | [Article]() | [WebServer]() | [Code Package]() | [Executable Program]()


Introduction
---------------
Nepre is a scoring function for calculating the protein's potential energy. And it will help to predict the 3D structure of a protein.
Nepre is designed to be **efficient**, **flexible** and **protable** with **two** typical algorithm included:
* Calculate the potential energy using different cutoff (Nepre-F).
* Calculate the potential energy using the radius of amino acid (Nepre-R). 

#### Nepre-F
Users can use the cutoff between 4 angstrom and 10 angstrom to calculate the potential energy.
#### Nepre-R  
30,000 high-resolution protein are used to get statistical result of radius and fitted by the gaussian distribution function. Then Nepre will calculate the potential energy according to the gaussian mean data of radius of each kind of amino acid.

CopyRight
-------------
**Nepre** is created by **LiuLab** of **Beijing Computation Science Research Center(CSRC)**.

Usage
-------------
Nepre is built completely on **python2.7**. Some basic
modules are required for using Nepre.  
* numpy
* tqmd

Contact Us
-------------
**Email**: nepre2018@163.com  
**Address**: Building 9, East Zone, ZPark II, No.10 East Xibeiwang Road, Haidian District, Beijing 100193, China.

Others
-------------
<div>Logo通过<a href="https://www.designevo.com/cn/" title="免费在线logo制作软件">DesignEvo</a>设计制作</div>