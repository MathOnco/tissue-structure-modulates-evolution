package Model;

import Framework.GridsAndAgents.AgentGrid2D;
import Framework.GridsAndAgents.AgentSQ2Dunstackable;
import Framework.Gui.*;
import Framework.Tools.FileIO;
import Framework.Rand;
import static Framework.Util.*;
import java.util.ArrayList;


class Parameters2D {
    public double mu = 1e-8;
    public double sp = 1e-3;
    public double sd = 0.1;
    public double Tp = 5e6;
    public double Td = 700;
    public double birth_rate = 0.5;
    public double death_rate = 0.5;
    public int r0 = 10;
    int sideLen = 200;
    int delete_thresh = 75; // ignore clone sizes smaller than this in Muller plots
}


class Cell2D extends AgentSQ2Dunstackable<PD2D> {
    int kd;
    int kp;
    int progenyID;
    int parentID;

    Cell2D Init(int kp0, int kd0, int progenyID0, int parentID0){
        kp = kp0;
        kd = kd0;
        progenyID = progenyID0;
        parentID = parentID0;
        G.driver_status[0] = 1;
        return this;
    }

    Cell2D Mutate(){
        boolean mutated = false;
        boolean driver_mutated = false;

        // driver mutation
        if((G.rn.Double() < ( G.Td * G.mu))) {
            kd++;
            if (kd > G.KdMAX) { G.KdMAX++; }
            mutated = true;
            driver_mutated = true;
        }

        // passenger mutation
        if((G.rn.Double() <(G.Tp * G.mu))) {
            kp++;
            if (kp > G.KpMAX) { G.KpMAX++; }
            mutated = true;
        }

        if (mutated) {
            parentID = progenyID;
            progenyID = G.progenyNextID;
            G.progenyToParentIDs[progenyID] = parentID;
            G.driver_status[progenyID] = (driver_mutated) ? kd : kd;
            G.progenyNextID++;
        }

        return this;
    }

    Cell2D Divide(){
        int nDivOptions = G.MapEmptyHood(G.neighborhood,Xsq(),Ysq());
        if(nDivOptions==0){ return null; }
        int nextAgentID = G.neighborhood[G.rn.Int(nDivOptions)];
        double first = (nextAgentID/G.sideLen) - (G.sideLen/2);
        double second = (nextAgentID%G.sideLen) - (G.sideLen/2);

        if ((first*first + second*second) <= (G.confRadius*G.confRadius)) {
            return G.NewAgentSQ(nextAgentID).Init(this.kp,this.kd,this.progenyID,this.parentID).Mutate();
        } else {
            return this;
        }
    }

    void Step(){

        // Passengers lower birth rate; Drivers raise birth rate
        if(G.rn.Double()<(Math.pow(1.0+G.sd,(double)kd)/Math.pow(1.0+G.sp,(double)kp)*G.BIRTH_RATE)){
            Divide();
        }

        // constant death rate
        if(G.rn.Double()<(G.DEATH_RATE )){
            Dispose();
            return;
        }
    }
}

public class PD2D extends AgentGrid2D<Cell2D> {

    // arrays used to store parentIDs for building Muller plots (and associated driver numbers)
    public int[] progenyToParentIDs = new int[7000000];
    public int[] driver_status = new int[7000000];

    // parameters
    public double mu = 1e-8;
    public double sp = 1e-3;
    public double sd = 0.1;
    public double Tp = 5e6;
    public double Td = 700;
    public double DEATH_RATE = 0.5;
    public double BIRTH_RATE = 0.5;
    public int r0 = 10;
    public int confRadius = 0;
    public int sideLen = 0;

    // tracking variables
    public int KpMAX = 0;
    public int KdMAX = 1;
    public int progenyNextID = 2; // assumes model is initialized with all cell's progenyID = 1


    // neighborhoods
    int[]neighborhood=MooreHood(false);
    Rand rn=new Rand(1);

    // constructor from parameters
    PD2D(Parameters2D p){
        //this(sideLen,r0);
        super(p.sideLen,p.sideLen, Cell2D.class,false,false);
        for (int x = 0; x < r0; x++) {
            for (int y = 0; y < r0; y++) {
                NewAgentSQ(x + p.sideLen/2 - r0 / 2,y + p.sideLen/2 - r0 / 2).Init(0,1,1,0);
            }
        }

        this.sideLen = p.sideLen;
        this.r0 = p.r0;
        this.mu = p.mu;
        this.sp = p.sp;
        this.sd = p.sd;
        this.Tp = p.Tp;
        this.Td = p.Td;
        this.BIRTH_RATE = p.birth_rate;
        this.DEATH_RATE = p.death_rate;
        this.confRadius = p.sideLen/2;

        // iterators / indicators
        this.KpMAX = 0;
        this.KdMAX = 1;

    }

    // third constructor (evenly distributed)
    PD2D(Parameters2D p, int regionSideLen, int density_modifier){
        super(regionSideLen,regionSideLen, Cell2D.class,false,false);

        for (int i = 0; i < regionSideLen*regionSideLen; i++) {
            if (i % density_modifier == 0) { NewAgentSQ(i).Init(0,1,1,0); }
        }

        this.sideLen = regionSideLen;
        this.r0 = p.r0;
        this.mu = p.mu;
        this.sp = p.sp;
        this.sd = p.sd;
        this.Tp = p.Tp;
        this.Td = p.Td;
        this.BIRTH_RATE = p.birth_rate;
        this.DEATH_RATE = p.death_rate;
        this.confRadius = 10000000; // "infinite" radius


        // iterators / indicators
        this.KpMAX = 0;
        this.KdMAX = 1;

    }

    // fourth constructor (used to initialize square models with dispersal)
    PD2D(Parameters2D p, int regionSideLen, int r0, boolean mullerBool){
        super(regionSideLen,regionSideLen, Cell2D.class,false,false);

        for (int x = 0; x < r0; x++) {
            for (int y = 0; y < r0; y++) {
                NewAgentSQ(x + regionSideLen/2 - r0 / 2,y + regionSideLen/2 - r0 / 2).Init(0,1,1,0);
            }
        }


        this.sideLen = regionSideLen;
        this.r0 = p.r0;
        this.mu = p.mu;
        this.sp = p.sp;
        this.sd = p.sd;
        this.Tp = p.Tp;
        this.Td = p.Td;
        this.BIRTH_RATE = p.birth_rate;
        this.DEATH_RATE = p.death_rate;
        this.confRadius = 10000000; // radius = "infinitely" large


        // iterators / indicators
        this.KpMAX = 0;
        this.KdMAX = 1;

    }


    // Step function for PD model ("steps" all cells through birth/death/mutation)
    void OriginalStep(){
        for (Cell2D c:this) {
            c.Step();
        }
        CleanShuffle(rn);
    }

    public static void SingleSim(Parameters2D p, int totalTime, int modifier, int sim, boolean saveTimelinesBoolean, String foldername) {

        PD2D model = new PD2D(p);



        // VISUALIZE
        UIWindow win = new UIWindow("Window", true);
        UIGrid Vis = new UIGrid(model.sideLen, model.sideLen, 1);
        TickTimer tt = new TickTimer();
        win.AddCol(0, Vis);
        win.RunGui();

        // INITIALIZE ARRAY
        // 1. time 2. pop 3. pass diversity 4. driver diversity
        // 5. n passengers 6. n drivers 7. max passengers 8. max drivers
        double[][] everything = new double[8][totalTime/modifier];
        for(int iii=0; iii<8; iii++) {
            for(int jjj=0; jjj<(totalTime/modifier); jjj++) {

                if (iii == 0) {
                    everything[iii][jjj] = (double)modifier * jjj;
                } else {
                    everything[iii][jjj] = 0.0;
                }
            }
        }

        int expectedProgeny = model.progenyToParentIDs.length;
        int[][] mullerDrivers = new int[20][totalTime / modifier];      // track the first 20 drivers (divisions)
        int[][] mullerGenetic = new int[expectedProgeny][totalTime / modifier]; // track the first [expectedProgeny] genetic clones
        for (int jjj = 0; jjj < (totalTime / modifier); jjj++) {
            for (int iii = 0; iii < 20; iii++) {
                mullerDrivers[iii][jjj] = 0;
            }
            for (int iii = 0; iii < expectedProgeny; iii++) {
                mullerGenetic[iii][jjj] = 0;
            }
        }

        // new gif
        String baseFilename =  foldername + "confRadius" + Integer.toString((int) (model.confRadius)) + "br" + Integer.toString((int) (model.BIRTH_RATE*100)) + "dr" + Integer.toString((int) (model.DEATH_RATE*100)) + "sp" + Integer.toString((int) (model.sp * 100000)) + "_sd" + Integer.toString((int) (model.sd * 100000)) + "_Tp" + Integer.toString((int) (model.Tp )) + "_Td" + Integer.toString((int) (model.Td ));
        FileIO everythingOutputFile =  new FileIO((baseFilename + ".csv"),( (sim == 0) ? "w" : "a" ));
        String gifFilename = baseFilename + Integer.toString(sim) + ".gif";
        GifMaker myGif = new GifMaker(gifFilename, 100,true);

        int j = 0;
        boolean keepStepping = true;
        for (int i = 0; i < totalTime; i++) {

            if (keepStepping) { model.OriginalStep(); }
            if (model.Pop() == 0) { keepStepping = false; }

            if (i % 50 == 0) {
                DrawCells(model, Vis);
                tt.TickPause(0);
                myGif.AddFrame(Vis);
                System.out.println("time: " + i + "pop: " + model.Pop());
            }

            // WRITE OUT VALUES
            if (i % modifier == 0) {

                if (saveTimelinesBoolean) {
                    everything[1][j] = model.Pop(); // total pop size
                    everything[2][j] = GetDiversity(model, 1); // pass div
                    everything[3][j] = GetDiversity(model, 0); // driv div
                    everything[4][j] = CountMutators(model, 1); // pass num
                    everything[5][j] = CountMutators(model, 0); // drivers num
                    everything[6][j] = model.KpMAX; // max pass
                    everything[7][j] = model.KdMAX; // max driv

                    // add all the Muller information
                    for (int k = 0; k < (model.sideLen*model.sideLen); k++) {
                        Cell2D c = model.GetAgent(k);
                        if (c != null) {
                            mullerDrivers[c.kd][j] = mullerDrivers[c.kd][j] + 1;
                            mullerGenetic[c.progenyID][j] = mullerGenetic[c.progenyID][j] + 1;
                        }
                    }
                }

                j++;
            }
        }

        if (saveTimelinesBoolean) {

            System.out.println("Building everything to output...");

            // OUTPUT "EVERYTHING" ARRAYS TO FILE
            BuildEverythingSingle(everythingOutputFile, everything, totalTime, modifier);

            // output driver parents (this is not calculated directly, but inferred from my logical tree)
            FileIO parentsDriverOutputFile = new FileIO((foldername + "driverParents.csv"), "w");
            WriteDriverParents(parentsDriverOutputFile);

            // OUTPUT DRIVER MULLER PLOTS:
            FileIO mullerDivOutputFile = new FileIO((foldername + "driverClones.csv"), "w");
            WriteDriverClones(mullerDivOutputFile,mullerDrivers,totalTime,modifier);

            System.out.println("Reducing phylogenies (this may take a minute)...");

            ////// genetic clones of REDUCED, single
            FileIO bigParentsReduced = new FileIO((foldername + "parents.csv"), "w");
            FileIO bigGeneticMullerReduced = new FileIO((foldername + "clones.csv"), "w");
            FileIO bigDriverStatusReduced = new FileIO((foldername + "driverStatus.csv"), "w");  //ryan
            ReduceParents(bigParentsReduced, bigGeneticMullerReduced, bigDriverStatusReduced, mullerGenetic, model.progenyToParentIDs, model.driver_status, expectedProgeny, totalTime, modifier, p.delete_thresh);

            System.out.println("Simulation finished...");
            System.out.println("Okay to close...");

        }

        win.Close();
    }

    /*
        Heatmaps()

        - iterates through mutation rate (mu) and passenger fitness penalty (sp)
        - outputs a matrix of number of number of non-extinct runs after totalTime
        - (such that N > 0 @ t = totalTime)
        - No visualizations, just data output at the end.
     */

    public static void Heatmaps(Parameters2D p, int nSims, int totalTime, String foldername) {

        int confRadius = p.sideLen / 2;
        String baseFilename =  foldername + "confRadius" + Integer.toString((int) (confRadius)) + "br" + Integer.toString((int) (p.birth_rate*100)) + "dr" + Integer.toString((int) (p.death_rate*100)) + "_sd" + Integer.toString((int) (p.sd * 100000)) + "_Tp" + Integer.toString((int) (p.Tp )) + "_Td" + Integer.toString((int) (p.Td ));

        int arrayLength = 0; int arrayLength2 = 0;
        for (double sp_power = -4.0; sp_power <= -1.0; sp_power+= 0.25) { arrayLength++; }
        for (double mu_power = -9.0; mu_power <= -6.0; mu_power+= 0.25) { arrayLength2++; }
        System.out.println("Planning to simulate " + arrayLength + " x " + arrayLength2 + " = " + arrayLength*arrayLength2 + " total simulations... ");


        FileIO successOutputFile =  new FileIO((baseFilename + ".csv"),"w");

        double mu_power;
        for (mu_power = -9.0; mu_power <= -6.0; mu_power+= 0.25) {

            successOutputFile =  new FileIO((baseFilename + ".csv"),"a");
            int[][] numSuccesses = new int[nSims][arrayLength]; // length of columns, one for each sp
            int[] numSuccessesSquish = new int[arrayLength]; // length of columns, one for each sp
            double mu_p = mu_power;

            System.out.println("Initiating multiple threads for mu = " + Math.pow(10, mu_p) + " for full range of sp ... ");
            MultiThread(nSims,8,(iThread)->{
                int arrayIter = 0;
                for (double sp_power = -4.0; sp_power <= -1.0; sp_power+= 0.25) {

                    double sp = Math.pow(10, sp_power);
                    double mu = Math.pow(10, mu_p);

                    PD2D model = new PD2D(p);

                    boolean keepStepping = true;
                    for (int i = 0; i < totalTime; i++) {
                        if (keepStepping) { model.OriginalStep(); }
                        if (model.Pop() == 0) { keepStepping = false; }
                    }

                    numSuccesses[iThread][arrayIter] = (model.Pop() > 0) ? 1 : 0;
                    arrayIter++;

                }
            });

            for (int i = 0; i < arrayLength; i++) {
                for (int j = 0; j < nSims; j++) {
                    numSuccessesSquish[i] += numSuccesses[j][i];
                }
            }

            successOutputFile.Write(ArrToString(numSuccessesSquish,",")+"\n");
            successOutputFile.Close();
        }
    }


    /*
        SegregatedRegions()

        - two simulations:
        - 1. Single region with a homogeneous initial condition (50% density of kd = 1 cells)
        - 2. Multiple regions: segregated into totalDivs x totalDivs regions
        - Note: this function doesn't save Muller plot information since there isn't a single founding cell
     */

    public static void SegregatedRegions(Parameters2D p, int totalTime, int modifier, int numberOfRowsAndColsOfRegions, String foldername) {

        // region side length
        int regionSideLen = p.sideLen/numberOfRowsAndColsOfRegions; // usually 310 (this is the square length, i think)
        int scaleFactor = 1;
        int density_modifier = 3;

        // Save stuff
        String baseFilename =  foldername + "sidebyside_" + "br" + Integer.toString((int) (p.birth_rate*100)) + "dr" + Integer.toString((int) (p.death_rate*100)) + "sp" + Integer.toString((int) (p.sp * 100000)) + "_sd" + Integer.toString((int) (p.sd * 100000)) + "_Tp" + Integer.toString((int) (p.Tp )) + "_Td" + Integer.toString((int) (p.Td ));

        for (int sim = 0; sim < 5; sim++) {

            // all the little models
            ArrayList<PD2D> models = new ArrayList<PD2D>();
            for (int row = 0; row < numberOfRowsAndColsOfRegions; row++) {
                for (int col = 0; col < numberOfRowsAndColsOfRegions; col++) {

                    PD2D nextModel = new PD2D(p, regionSideLen, density_modifier);
                    models.add(nextModel);
                }
            }

            // single region model
            PD2D bigModel = new PD2D(p, p.sideLen, density_modifier);

            // OUTPUT FILE
            FileIO everythingBigOutputFile = new FileIO((baseFilename + "_" + Integer.toString(numberOfRowsAndColsOfRegions) + "DivsBig.csv"), ((sim == 0) ? "w" : "a"));
            FileIO everythingDivsOutputFile = new FileIO((baseFilename + "_" + Integer.toString(numberOfRowsAndColsOfRegions) + "Divs.csv"), ((sim == 0) ? "w" : "a"));


            // INITIALIZE BIG ARRAY
            // 1. time 2. pop 3. pass diversity 4. driver diversity
            // 5. n passengers 6. n drivers 7. max passengers 8. max drivers
            double[][] everythingBig = new double[8][totalTime / modifier];
            for (int iii = 0; iii < 8; iii++) {
                for (int jjj = 0; jjj < (totalTime / modifier); jjj++) {

                    if (iii == 0) {
                        everythingBig[iii][jjj] = (double) modifier * iii;
                    } else {
                        everythingBig[iii][jjj] = 0.0;
                    }
                }
            }

            // INITIALIZE DIVS ARRAY
            // 1. time 2. pop 3. pass diversity 4. driver diversity
            // 5. n passengers 6. n drivers 7. max passengers 8. max drivers
            double[][] everythingDivs = new double[8][totalTime / modifier];
            for (int iii = 0; iii < 8; iii++) {
                for (int jjj = 0; jjj < (totalTime / modifier); jjj++) {

                    if (iii == 0) {
                        everythingDivs[iii][jjj] = (double) modifier * iii;
                    } else {
                        everythingDivs[iii][jjj] = 0.0;
                    }
                }
            }


            // VISUALIZE
            UIWindow win = new UIWindow("Window", true);
            UIGrid Vis = new UIGrid(regionSideLen * numberOfRowsAndColsOfRegions, regionSideLen * numberOfRowsAndColsOfRegions, scaleFactor);
            win.AddCol(0, Vis);
            win.RunGui();

            // the other one (big window)
            UIWindow BigWin = new UIWindow("Big Window", true);
            UIGrid BigViz = new UIGrid(regionSideLen * numberOfRowsAndColsOfRegions, regionSideLen * numberOfRowsAndColsOfRegions, scaleFactor);
            BigWin.AddCol(0, BigViz);
            BigWin.RunGui();

            // new gif
            GifMaker smallGif = new GifMaker(baseFilename + "_" + Integer.toString(numberOfRowsAndColsOfRegions) + "Divs.gif", 100, true);
            GifMaker bigGif = new GifMaker(baseFilename + "_" + Integer.toString(numberOfRowsAndColsOfRegions) + "DivsBig.gif", 100, true);

            int j = 0;

            boolean keepStepping = true;
            for (int i = 0; i < totalTime; i++) {

                // step all models
                if (keepStepping) {
                    for (int k = 0; k < models.size(); k++) {
                        models.get(k).OriginalStep();
                    }
                    bigModel.OriginalStep();
                }


                // WRITE OUT VALUES
                if (i % modifier == 0) {
                    //System.out.println(sim + " : " + model.confRadius + " : " +  i);
                    System.out.println(numberOfRowsAndColsOfRegions + " divs & time: " + i);

                    for (int rowi = 0; rowi < numberOfRowsAndColsOfRegions; rowi++) {
                        for (int coli = 0; coli < numberOfRowsAndColsOfRegions; coli++) {
                            DrawOffsetCells(models.get(rowi + numberOfRowsAndColsOfRegions * coli), Vis, rowi, coli);
                        }
                    }

                    DrawCells(bigModel, BigViz);
                    smallGif.AddFrame(Vis);
                    bigGif.AddFrame(BigViz);


                    // SAVE everythingBig
                    everythingBig[1][j] = bigModel.Pop(); // total pop size
                    everythingBig[2][j] = GetDiversity(bigModel, 1); // pass div
                    everythingBig[3][j] = GetDiversity(bigModel, 0); // driv div
                    everythingBig[4][j] = CountMutators(bigModel, 1); // pass num
                    everythingBig[5][j] = CountMutators(bigModel, 0); // drivers num
                    everythingBig[6][j] = bigModel.KpMAX; // max pass
                    everythingBig[7][j] = bigModel.KdMAX; // max driv

                    // SAVE everythingDivs NEED TO CHANGE
                    everythingDivs[1][j] = GetRegionPop(models); // total pop size
                    everythingDivs[2][j] = GetRegionDiversity(models, 1); // pass div
                    everythingDivs[3][j] = GetRegionDiversity(models, 0); // driv div
                    everythingDivs[4][j] = CountRegionMutators(models, 1); // pass num
                    everythingDivs[5][j] = CountRegionMutators(models, 0); // drivers num
                    everythingDivs[6][j] = GetRegionKMax(models, 0); // max pass
                    everythingDivs[7][j] = GetRegionKMax(models, 1); // max driv

                    j++;
                }
            }

            // OUTPUT TO FILE
            for (int jj = 0; jj < 8; jj++) {

                StringBuilder sb = new StringBuilder();
                for (int ii = 0; ii < (totalTime / modifier) - 1; ii++) {
                    sb.append(everythingBig[jj][ii] + ",");
                }
                sb.append(everythingBig[jj][(totalTime / modifier) - 1] + "\n");
                everythingBigOutputFile.Write(sb.toString());
            }
            everythingBigOutputFile.Close();

            // OUTPUT DIVS TO FILE
            for (int jj = 0; jj < 8; jj++) {

                StringBuilder sb = new StringBuilder();
                for (int ii = 0; ii < (totalTime / modifier) - 1; ii++) {
                    sb.append(everythingDivs[jj][ii] + ",");
                }
                sb.append(everythingDivs[jj][(totalTime / modifier) - 1] + "\n");
                everythingDivsOutputFile.Write(sb.toString());
            }
            everythingDivsOutputFile.Close();


            // CLOSE GIFS, WINDOWS
            smallGif.Close();
            bigGif.Close();

            win.Close();
            BigWin.Close();

        }

        return;
    }


    public static void SegregatedRegionsWithDispersal(Parameters2D p, int totalTime, int modifier, int numberOfRowsAndColsOfRegions, String foldername, double dispersalRate) {

        int regionSideLen = p.sideLen/numberOfRowsAndColsOfRegions;
        int scaleFactor = 1;

        // all the cute little baby models in a grid
        ArrayList<PD2D> models = new ArrayList<PD2D>();
        int expectedProgeny = 0;
        for (int row = 0; row < numberOfRowsAndColsOfRegions; row++) {
            for (int col = 0; col < numberOfRowsAndColsOfRegions; col++) {

                PD2D nextModel = new PD2D(p, regionSideLen, ((row == (int)(( numberOfRowsAndColsOfRegions - 1 ) / 2)) && (col == (int)(( numberOfRowsAndColsOfRegions - 1 ) / 2))) ? p.r0 : 0, true);
                expectedProgeny = nextModel.progenyToParentIDs.length;
                models.add(nextModel);
            }
        }

        // OUTPUT FILES
        String baseFilename =  foldername + "everything_" + "br" + Integer.toString((int) (p.birth_rate*100)) + "dr" + Integer.toString((int) (p.death_rate*100)) + "sp" + Integer.toString((int) (p.sp * 100000)) + "_sd" + Integer.toString((int) (p.sd * 100000)) + "_Tp" + Integer.toString((int) (p.Tp )) + "_Td" + Integer.toString((int) (p.Td ));
        FileIO everythingDivsOutputFile = new FileIO((baseFilename + "_" + Integer.toString(numberOfRowsAndColsOfRegions) + "Divs.csv"), "w");

        // INITIALIZE BIG & DIVS ARRAY
        // 1. time 2. pop 3. pass diversity 4. driver diversity 5. n passengers 6. n drivers 7. max passengers 8. max drivers
        double[][] everythingDivs = new double[8][totalTime / modifier];
        int[][] mullerDriversDivs = new int[20][totalTime / modifier];      // track the first 20 drivers (divisions)
        int[][] mullerGeneticDivs = new int[expectedProgeny][totalTime / modifier]; // track the first [expectedProgeny] genetic clones

        // initialize all arrays
        for (int jjj = 0; jjj < (totalTime / modifier); jjj++) {
            for (int iii = 0; iii < 8; iii++) {
                everythingDivs[iii][jjj] = (iii == 0) ? (double) modifier * iii : 0.0;
            }
        }
        for (int jjj = 0; jjj < (totalTime / modifier); jjj++) {
            for (int iii = 0; iii < 20; iii++) {
                mullerDriversDivs[iii][jjj] = 0;
            }
            for (int iii = 0; iii < expectedProgeny; iii++) {
                mullerGeneticDivs[iii][jjj] = 0;
            }
        }

        // VISUALIZE
        UIWindow win; UIWindow BigWin;  GifMaker smallGif; GifMaker bigGif; UIGrid Vis; UIGrid BigViz;
        win = new UIWindow("Window", true);
        smallGif = new GifMaker(foldername + "gifDiv.gif", 100, true);
        Vis = new UIGrid(regionSideLen * numberOfRowsAndColsOfRegions, regionSideLen * numberOfRowsAndColsOfRegions, scaleFactor);
        win.AddCol(0, Vis);
        win.RunGui();

        // re-use for finding empty indices (for dispersal between spatial divisions / bins)
        int[] emptyIndices = new int[p.sideLen*p.sideLen];
        int j = 0;
        int globalProgenyNextID = 1; // for smalls/divs models

        for (int i = 0; i < totalTime; i++) {
            if (i % 50 == 0) {
                for (int rowi = 0; rowi < numberOfRowsAndColsOfRegions; rowi++) {
                    for (int coli = 0; coli < numberOfRowsAndColsOfRegions; coli++) {
                        DrawOffsetCells(models.get(rowi + numberOfRowsAndColsOfRegions * coli), Vis, rowi, coli);
                    }
                }

                smallGif.AddFrame(Vis);
            }

            // WRITE OUT VALUES
            if (i % modifier == 0) {
//                System.out.println(numberOfRowsAndColsOfRegions + " divs & time: " + i);

                // SAVE everythingDivs
                everythingDivs[1][j] = GetRegionPop(models); // total pop size
                everythingDivs[2][j] = GetRegionDiversity(models, 1); // pass div
                everythingDivs[3][j] = GetRegionDiversity(models, 0); // driv div
                everythingDivs[4][j] = CountRegionMutators(models, 1); // pass num
                everythingDivs[5][j] = CountRegionMutators(models, 0); // drivers num
                everythingDivs[6][j] = GetRegionKMax(models, 0); // max pass
                everythingDivs[7][j] = GetRegionKMax(models, 1); // max driv


                // small models (with divisions) save muller information:
                for (int rowi = 0; rowi < numberOfRowsAndColsOfRegions; rowi++) {
                    for (int coli = 0; coli < numberOfRowsAndColsOfRegions; coli++) {
                        PD2D nthModel = models.get(rowi + numberOfRowsAndColsOfRegions * coli);
                        for (int k = 0; k < (regionSideLen*regionSideLen); k++) {
                            Cell2D c = nthModel.GetAgent(k);
                            if (c != null) {
                                mullerDriversDivs[c.kd][j] = mullerDriversDivs[c.kd][j] + 1;
                                mullerGeneticDivs[c.progenyID][j] = mullerGeneticDivs[c.progenyID][j] + 1;
                            }
                        }
                    }
                }
                j++;
            }

            // check for dispersal / circulation events
            for (int rowk = 0; rowk < numberOfRowsAndColsOfRegions; rowk++) {
                for (int colk = 0; colk < numberOfRowsAndColsOfRegions; colk++) {
                    PD2D focusModel = models.get(rowk + numberOfRowsAndColsOfRegions * colk);
                    double effCircRate = dispersalRate*(double)(focusModel.Pop()/(double)(regionSideLen*regionSideLen));
                    if ( ( focusModel.rn.Double() < effCircRate ) && (focusModel.Pop() > 0 )){

                        // which other model is it going to?
                        int rowIter = focusModel.rn.Int(3) - 1;
                        int colIter = focusModel.rn.Int(3) - 1;
                        int rowNext = rowIter + rowk;
                        int colNext = colIter + colk;
                        if (rowNext > (numberOfRowsAndColsOfRegions-1)) { rowNext = 0; }
                        if (rowNext < 0) { rowNext = numberOfRowsAndColsOfRegions - 1; }
                        if (colNext > (numberOfRowsAndColsOfRegions-1)) { colNext = 0; }
                        if (colNext < 0) { colNext = numberOfRowsAndColsOfRegions - 1; }

                        // random leaving agent
                        Cell2D leavingAgent = focusModel.RandomAgent(focusModel.rn);

                        // random starting point
                        PD2D circModel = models.get(rowNext + numberOfRowsAndColsOfRegions * colNext);

                        int lastAgent = 0;
                        for (int empty = 0; empty < regionSideLen*regionSideLen; empty++) {
                            if (circModel.GetAgent(empty) == null) {
                                emptyIndices[lastAgent] = empty;
                                lastAgent++;
                            }
                        }

                        int searchingIndex = emptyIndices[focusModel.rn.Int(lastAgent-1)];

                        // focusModel.leavingAgent -> circModel(searchingIndex)
                        circModel.NewAgentSQ(searchingIndex).Init(leavingAgent.kp,leavingAgent.kd, leavingAgent.progenyID, leavingAgent.parentID);
                        focusModel.GetAgent(leavingAgent.Isq()).Dispose();
                    }
                }
            }

            // step all models
            for (int k = 0; k < models.size(); k++) {
                models.get(k).OriginalStep();

                // check if all models have synchronized progenyNextID here!!
                if (models.get(k).progenyNextID > globalProgenyNextID) {
                    // update all others
                    globalProgenyNextID = models.get(k).progenyNextID;
                    for (int kk = 0; kk < models.size(); kk++) {
                        models.get(kk).progenyNextID = globalProgenyNextID;
                    }
                }

            }
        }


        System.out.println("Building everything to output...");

        // OUTPUT "EVERYTHING" ARRAYS TO FILE
        BuildEverythingSingle(everythingDivsOutputFile, everythingDivs, totalTime, modifier);


        // output driver parents (this is not calculated directly, but inferred from my logical tree)
        FileIO parentsDriverDivOutputFile = new FileIO((foldername + "driverParentsDiv.csv"), "w");
        //WriteDriverParents(parentsDriverBigOutputFile);
        WriteDriverParents(parentsDriverDivOutputFile);

        // OUTPUT DRIVER MULLER PLOTS:
        FileIO mullerDivOutputFile = new FileIO((foldername + "driverClonesDiv.csv"), "w");
        //WriteDriverClones(mullerBigOutputFile,mullerDriversBig,totalTime,modifier);
        WriteDriverClones(mullerDivOutputFile,mullerDriversDivs,totalTime,modifier);

        // combine div parents and driver status into a global vector
        int globalProgenyToParentIDs[] = new int[expectedProgeny];
        int globalDriverStatus[] = new int[expectedProgeny]; // ryan
        for (int i = 1; i < expectedProgeny; i++) {
            for (int kk = 0; kk < models.size(); kk++) {
                // change global vector only if this model is non-zero
                globalProgenyToParentIDs[i] = (models.get(kk).progenyToParentIDs[i] != 0) ? models.get(kk).progenyToParentIDs[i] : globalProgenyToParentIDs[i];
                globalDriverStatus[i] = (models.get(kk).driver_status[i] != 0) ? models.get(kk).driver_status[i] : globalDriverStatus[i];
            }
        }

        // now we have two parent vectors: globalProgenyToParentIDs[ii] and bigModel.progenyToParentIDs[ii]
        // we also have two driver status vectors: globalDriverStatus[ii] and bigModel.driver_status[ii]

        System.out.println("Reducing phylogenies (this may take a minute)...");
        FileIO divsParentsReduced = new FileIO((foldername + "parentsDiv.csv"), "w");
        FileIO divsGeneticMullerReduced = new FileIO((foldername + "clonesDiv.csv"), "w");
        FileIO divsDriverStatusReduced = new FileIO((foldername + "driverStatusDiv.csv"), "w"); //ryan
        ReduceParents(divsParentsReduced,divsGeneticMullerReduced,divsDriverStatusReduced,mullerGeneticDivs, globalProgenyToParentIDs, globalDriverStatus, expectedProgeny, totalTime, modifier, p.delete_thresh);

        System.out.println("Simulation finished...");
        System.out.println("Okay to close...");

        // CLOSE GIFS, WINDOWS
        smallGif.Close();
        win.Close();

        return;
    }

    public static void main(String[] args) {

        /*

            SINGLE SIMULATION, constrained to a circular domain

                - uncomment out the following section to run a single sim
                - outputs basic information ("everything" array) and GIFs and Muller plot information

         */

        int modifier = 50; // when to save data
        int totalTime = 1000; // total simulation time
        String foldername = "data-output/";
        Parameters2D parameters = new Parameters2D();
        parameters.delete_thresh = 20;
        SingleSim(parameters,totalTime,modifier,0,true,foldername);



        /*

            HEATMAP of tumor non-extinction events

                - uncomment out the following section to get a heatmap of tumor extinction events
                - runs simulations for given parameters and a range of sp and mu values
                - extinction events are when N = 0 after totalTime time steps
                - outputs a matrix of number of simulations for which tumor does not go extinct

         */

//        int nSims = 8; // when to save data
//        int totalTime = 1000; // total simulation time
//        String foldername = "data-output/";
//        Parameters2D parameters = new Parameters2D();
//        Heatmaps(parameters, nSims,totalTime, foldername);



        /*

            SEGREGATED REGIONS

                - uncomment out the following section to run two simulations:
                    1. a single square domain with homogeneous initial condition (33% density)
                    2. same size domain, segregated into numberOfRowsAndColsOfRegions^2 regions
                - each segregated region does not interact with each other region.
                - outputs basic statistics ("everything" arrays) and GIFs

         */

//        int modifier = 50; // when to save data
//        int totalTime = 1000; // total simulation time
//        String foldername = "data-output/";
//        Parameters2D parameters = new Parameters2D();
//        int numberOfRowsAndColsOfRegions = 5; // (i.e. 4 by 4 regions)
//        parameters.sideLen = 500;
//        SegregatedRegions(parameters,totalTime,modifier,numberOfRowsAndColsOfRegions, "data-output/");









        /*

            SEGREGATED REGIONS WITH DISPERSAL

                - uncomment out the following section to run two simulations:
                    1. a single square domain with homogeneous initial condition (33% density)
                    2. same size domain, segregated into numberOfRowsAndColsOfRegions^2 regions
                - each segregated region interacts with each other region through dispersal events
                    - dispersal events occur at a rate of dispersalRate*cellDensity
                - outputs basic statistics ("everything" arrays) and GIFs

         */

//        int modifier = 100; // when to save data
//        int totalTime = 1000; // total simulation time
//        String foldername = "data-output/";
//        Parameters2D parameters = new Parameters2D();
//        int numberOfRowsAndColsOfRegions = 5; // (i.e. 4 by 4 regions)
//        parameters.sideLen = 500;
//        double dispersalRate = 0.9;
//        SegregatedRegionsWithDispersal(parameters,totalTime,modifier,numberOfRowsAndColsOfRegions,foldername,dispersalRate);
//




        return;
    }

    public static void DrawCells(PD2D model, UIGrid visCells) {
        // color half by drivers and half by passengers
        for (int i = 0; i < visCells.length; i++) {
            Cell2D c=model.GetAgent(i);

            if(c==null){
                visCells.SetPix(i, RGB(1,1,1));
            } else{
                visCells.SetPix(i,()->{return CategorialColor((c.kd - 1) % 19 );});
            }
        }
    }

    public static void DrawOffsetCells(PD2D model, UIGrid visCells, int row, int col) {
        // color half by drivers and half by passengers
        for (int x = 0; x < model.sideLen; x++) {
            for (int y = 0; y < model.sideLen; y++) {

                Cell2D c = model.GetAgent(x,y);

                int shift_x = x + model.sideLen*col;
                int shift_y = y + model.sideLen*row;

                if (c== null) {
                    visCells.SetPix(shift_x,shift_y,RGB(1,1,1));
                } else {
                    visCells.SetPix(shift_x,shift_y,()->{return CategorialColor((c.kd - 1) % 19 );});
                }

            }
        }
    }

    public static double GetDiversity(PD2D model, int passBool) {
        if (passBool == 0) {
            // count drivers!
            int kVec[] = new int[model.KdMAX + 1];
            for (Cell2D c : model) {
                kVec[c.kd - 1]++;
            }

            // calculate diversity
            double sum = 0;
            for (int i = 0; i <= model.KdMAX; i++) {
                if (kVec[i] > 0) {
                    sum += (double)((double)kVec[i] / (double)model.Pop())*Math.log10((double)kVec[i] / (double)model.Pop());
                }
            }
            return Math.exp(- sum);

        } else {
            // count passengers!
            int kVec[] = new int[model.KpMAX + 1];
            for (Cell2D c : model) {
                kVec[c.kp]++;
            }

            // calculate diversity
            double sum = 0;
            for (int i = 0; i <= model.KpMAX; i++) {
                if (kVec[i] > 0) {
                    sum += (double)((double)kVec[i] / (double)model.Pop())*Math.log10((double)kVec[i] / (double)model.Pop());
                }
            }
            return Math.exp(- sum);

        }
    }

    public static double GetRegionDiversity(ArrayList<PD2D> models, int passBool) {

        int totalKpMax = 0;
        int totalKdMax = 0;
        int totalPop = 0;

        for (int k = 0; k < models.size(); k++) {
            PD2D thisModel = models.get(k);
            totalKpMax = (thisModel.KpMAX > totalKpMax) ? thisModel.KpMAX : totalKpMax;
            totalKdMax = (thisModel.KdMAX > totalKdMax) ? thisModel.KdMAX : totalKdMax;
            totalPop += thisModel.Pop();
        }


        if (passBool == 0) {
            // count drivers!
            int kVec[] = new int[totalKdMax + 1];

            for (int k = 0; k < models.size(); k++) {
                PD2D thisModel = models.get(k);
                for (Cell2D c : thisModel) {
                    kVec[c.kd - 1]++;
                }
            }

            // calculate diversity
            double sum = 0;
            for (int i = 0; i <= totalKdMax; i++) {
                if (kVec[i] > 0) {
                    sum += (double)((double)kVec[i] / (double)totalPop)*Math.log10((double)kVec[i] / (double)totalPop);
                }
            }
            return Math.exp(- sum);

        } else {
            // count passengers!
            int kVec[] = new int[totalKpMax + 1];


            for (int k = 0; k < models.size(); k++) {
                PD2D thisModel = models.get(k);
                for (Cell2D c : thisModel) {
                    kVec[c.kp]++;
                }
            }

            // calculate diversity
            double sum = 0;
            for (int i = 0; i <= totalKpMax; i++) {
                if (kVec[i] > 0) {
                    sum += (double)((double)kVec[i] / (double)totalPop)*Math.log10((double)kVec[i] / (double)totalPop);
                }
            }
            return Math.exp(- sum);
        }
    }

    public static double GetRegionKMax(ArrayList<PD2D> models, int passBool) {

        int totalKpMax = 0;
        int totalKdMax = 0;

        for (int k = 0; k < models.size(); k++) {
            PD2D thisModel = models.get(k);
            totalKpMax = (thisModel.KpMAX > totalKpMax) ? thisModel.KpMAX : totalKpMax;
            totalKdMax = (thisModel.KdMAX > totalKdMax) ? thisModel.KdMAX : totalKdMax;
        }

        return (passBool == 0) ? totalKpMax : totalKdMax;

    }

    public static double GetRegionPop(ArrayList<PD2D> models) {

        int totalPop = 0;

        for (int k = 0; k < models.size(); k++) {
            PD2D thisModel = models.get(k);
            totalPop += thisModel.Pop();
        }

        return totalPop;

    }

    public static int CountMutators(PD2D model, int passBool) {
        int sum = 0;
        if (passBool == 0) {
            // count drivers!

            for (Cell2D c : model) {
                sum += c.kd;
            }
            return sum;

        } else {
            // count passengers!
            for (Cell2D c : model) {
                sum += c.kp;
            }

            return sum;
        }


    }

    public static int CountRegionMutators(ArrayList<PD2D> models, int passBool) {
        int sum = 0;
        if (passBool == 0) {
            for (int k = 0; k < models.size(); k++) {
                for (Cell2D c : models.get(k)) {
                    sum += c.kd;
                }
            }

        } else {
            for (int k = 0; k < models.size(); k++) {
                for (Cell2D c : models.get(k)) {
                    sum += c.kp;
                }
            }
        }
        return sum;
    }

    public static void BuildEverythingSingle(FileIO everythingDivsOutputFile, double[][] everythingDivs, int totalTime, int modifier) {
        StringBuilder sb = new StringBuilder();

        for (int jj = 0; jj < 8; jj++) {
            sb = new StringBuilder();
            for (int ii = 0; ii < (totalTime / modifier) - 1; ii++) { sb.append(everythingDivs[jj][ii] + ","); }
            sb.append(everythingDivs[jj][(totalTime / modifier) - 1] + "\n");
            everythingDivsOutputFile.Write(sb.toString());
        }

        everythingDivsOutputFile.Close();
    }

    public static void WriteDriverParents(FileIO outputFile) {
        // write out parents beginning
        outputFile.Write(Integer.toString((int)0));


        for (int jj = 1; jj < 20; jj++) {
            outputFile.Write("," + Integer.toString((int)jj));
        }

        outputFile.Close();
    }

    public static void WriteDriverClones(FileIO outputFile, int[][] mullerDrivers, int totalTime, int modifier) {
        // first, output time:
        StringBuilder sb = new StringBuilder();

        for (int time = 0; time < (totalTime / modifier) - 1; time++) { sb.append(time*modifier + ","); }
        sb.append(modifier*((totalTime / modifier) - 1) + "\n");
        outputFile.Write(sb.toString());

        for (int jj = 0; jj < 20; jj++) {
            sb = new StringBuilder();
            for (int ii = 0; ii < (totalTime / modifier) - 1; ii++) {
                sb.append(mullerDrivers[jj][ii] + ",");
            }
            sb.append(mullerDrivers[jj][(totalTime / modifier) - 1] + "\n");
            outputFile.Write(sb.toString());
        }

        outputFile.Close();
    }


    /*

        ReduceParents()
            - removes clones that are dead and gone
            - reconnects parents together so that there is a continuous lineage when parent dies and child lives on
            - writes everything out

     */

    public static void ReduceParents(FileIO parentsReduced, FileIO geneticMullerReduced, FileIO driverStatusReduced, int[][] mullerGenetic, int[] progenyToParentIDs, int[] driver_status, int expectedProgeny, int totalTime, int modifier, int delete_thresh) {


        int rowmax = 0;
        ArrayList<Integer> reducedParents = new ArrayList<>();
        ArrayList<Integer> reducedProgeny = new ArrayList<>();
        ArrayList<Integer> reducedDriverStatus = new ArrayList<>();
        reducedParents.add(new Integer(-1)); // add 0th parent
        reducedProgeny.add(new Integer(0)); // add 0th progeny

        for (int jj = 1; jj < expectedProgeny; jj++) {
            rowmax = 0;
            for (int ii = 0; ii < (totalTime / modifier) - 1; ii++) {
                rowmax = (mullerGenetic[jj][ii] > rowmax) ? mullerGenetic[jj][ii]: rowmax;
            }

            if (rowmax > delete_thresh) {

                int mySupposedParent = progenyToParentIDs[jj];

                // check if my parent is in there:
                while (!reducedProgeny.contains(new Integer(mySupposedParent))) {

                    // parent is parent of my parent, lol -- is that parent in there?
                    mySupposedParent = progenyToParentIDs[progenyToParentIDs[mySupposedParent]];

                    // breaks if it is in there, eventually up the chain
                }

                // if my parent is in the reduced thing, add me.
                reducedParents.add(new Integer(reducedProgeny.indexOf(mySupposedParent)));
                reducedProgeny.add(new Integer(jj));
                reducedDriverStatus.add(new Integer(driver_status[jj]));

            }
        }


        // first, output time:
        StringBuilder sb = new StringBuilder();
        StringBuilder sb2 = new StringBuilder();
        StringBuilder sb3 = new StringBuilder(); //ryan
        for (int time = 0; time < (totalTime / modifier) - 1; time++) { sb.append(time*modifier + ","); }
        sb.append(modifier*((totalTime / modifier) - 1) + "\n");
        geneticMullerReduced.Write(sb.toString()); // file name


        // out gene muller big to file
        int ii_reduced = 0;
        for (int jj = 1; jj < expectedProgeny; jj++) {
            sb = new StringBuilder();

            // check if it's in reduced
            if (reducedProgeny.contains(new Integer(jj))) {
                for (int ii = 0; ii < (totalTime / modifier) - 1; ii++) {
                    sb.append(mullerGenetic[jj][ii] + ",");
                }
                sb.append(mullerGenetic[jj][(totalTime / modifier) - 1] + "\n");

                // output string, including newline \n
                geneticMullerReduced.Write(sb.toString());


                // first iteration no comma
                sb2.append(( (ii_reduced == 0) ? Integer.toString(reducedParents.get(ii_reduced) + 1) : "," + Integer.toString(reducedParents.get(ii_reduced) + 1) )); // add supposed parent, not actual

                // first one is driver, otherwise check driver status
                sb3.append(( (ii_reduced == 0) ? Integer.toString(1) : "," + Integer.toString(reducedDriverStatus.get(ii_reduced)) ));

                ii_reduced++;
            }


        }
        parentsReduced.Write(sb2.toString());
        driverStatusReduced.Write(sb3.toString());
        driverStatusReduced.Close();
        geneticMullerReduced.Close();
        parentsReduced.Close();

    }

}

