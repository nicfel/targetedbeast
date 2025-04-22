## TargetedBeast pacakge for [BEAST 2](https://beast2.org)


## Install the package

### Installing through BEAUti

TargetedBeast is a [BEAST2](http://beast2.org) package that requires BEAST 2 v2.7.7.
If you have not already done so, you can get BEAST 2 from [here](http://beast2.org).

To install TargetedBeast, it is easiest to start BEAUti (a program that is part of BEAST), and select the menu `File/Manage packages`. A package manager dialog pops up, that looks something like this:

![Package Manager](https://github.com/CompEvol/CCD/raw/master/doc/package_repos.png)

If the TargetedBeast package is listed, just click on it to select it, and hit the `Install/Upgrade` button.

If the TargetedBeast package is not listed, you may need to add a package repository by clicking the `Package repositories` button. A window pops up where you can click `Add URL` and add `https://raw.githubusercontent.com/CompEvol/CBAN/master/packages-extra-2.7.xml` in the entry. After clicking OK, the dialog should look something like this:

![Package Repositories](https://github.com/CompEvol/CCD/raw/master/doc/package_repos0.png)

Click OK and now TargetedBeast should be listed in the package manager (as in the first dialog above). Select and click Install/Upgrade to install.

### Install by hand

* Download the package from [here](https://github.com/nicfel/targetedbeast/releases/download/v0.8.0/TargetedBeast.v0.8.0.zip)
* Create TargetedBeast directory inside BEAST package directory
  * for Windows in Users\<YourName>\BEAST\2.X\TargetedBeast
  * for Mac in /Users/<YourName>\/Library/Application Support/BEAST/2.X/TargetedBeast
  * for Linux /home/<YourName>/.beast/2.X/TargetedBeast
  Here <YourName> is the username you use, and in “2.X” the X refers to the major version of BEAST, so 2.X=2.7 for version 2.7.7.
* Unzip the file `TargetedBeast.v0.0.8.zip` inside the TargetedBeast directory

## Build from code

* Get code for beast2, BeastFX and CCD repositories:
  * git clone https://github.com/CompEvol/beast2.git
  * git clone https://github.com/CompEvol/BeastFX.git
  * git clone https://github.com/nicfel/TargetedBeast.git
* Run `ant install` from the TargetedBeast directory
  


## Using TargetedBeast

## In BEAUti

* Set up your analysis as per usual
* Select the menu `View=>Show Operators` panel

![Show operators](https://github.com/nicfel/targetedbeast/raw/master/doc/show-operators.png)


* Select TargetedOperatorSchedulefrom the drop-down box under Operator schedule at the bottom of the screen.

![Select targeted schedule](https://github.com/nicfel/targetedbeast/raw/master/doc/select-targeted-schedule.png)

After saving, the operator schedule is replaced, and when running the XML in BEAST a message appears after the list of citations that looks something like this:

```
replacing YuleModelSubtreeSlide.t:dna with BactrianSubtreeSlide and RangeSlide
replacing YuleModelNarrow.t:dna with WeightBasedNodeRandomizer and HeightBasedNodeRandomizer
replacing YuleModelWide.t:dna with WeightedWideOperator and WeightedWideOperator
replacing YuleModelWilsonBalding.t:dna with TargetedWilsonBalding and TargetedWilsonBaldingusing edge lengths
replacing YuleModelBICEPSEpochTop.t:dna with IntervalScaleOperator
removing YuleModelBICEPSEpochAll.t:dna (replaced by IntervalScaleOperator)
removing YuleModelBICEPSTreeFlex.t:dna (replaced by IntervalScaleOperator)
```

## By hand

Open the XML in a text editor and search/replace `spec="OperatorSchedule"` with `spec="targetedbeast.operatorschedule.TargetedOperatorSchedule"`, save the file and run in BEAST.



