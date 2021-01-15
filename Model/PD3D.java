package Model;
import Framework.GridsAndAgents.AgentGrid3D;
import Framework.GridsAndAgents.AgentSQ3Dunstackable;
import Framework.GridsAndAgents.Grid3Dint;
import Framework.Gui.*;
import Framework.Tools.FileIO;
import Framework.Rand;
import static Framework.Util.*;
import java.util.ArrayList;
import Framework.Gui.OpenGL3DWindow;


class Parameters3D {
    public double mu = 1e-8;
    public double sp = 1e-3;
    public double sd = 0.1;
    public double Tp = 5e6;
    public double Td = 700;
    public double DEATH_RATE = 0.5;
    public double BIRTH_RATE = 0.5;
    public int mode = 1; // P-D birth-death mode
    public int r0 = 10;
    public final static int OPEN_DUCT=1,DUCT_WALL=0;

    // tracking variables
    public int KpMAX = 0;
    public int KdMAX = 1;
    public int progenyNextID = 2; // assume all are started with 1

    // more parameters
    public int delete_thresh = 100;
    public int totalTime = 1000;
    public int modifier = 20;
    public int initialCells = 100;
    public int kd0 = 1;

    public int sideX = 500;
    public int sideY = 500;
    public int sideZ = 100;

    public int progenyArrayLength = 3500000;

    public boolean INVASION = false;

}

class Cell3D extends AgentSQ3Dunstackable<PD3D> {
    int kd;
    int kp;
    int progenyID;
    int parentID;

    Cell3D Init(int kp0, int kd0, int progenyID0, int parentID0){
        kp = kp0;
        kd = kd0;
        progenyID = progenyID0;
        parentID = parentID0;
        return this;
    }
    Cell3D Mutate(){
        boolean mutated = false;
        boolean driver_mutated = false;
        // driver mutation
        if((G.rn.Double() < ( G.p.Td * G.p.mu))) {
            kd++;
            if (kd > G.p.KdMAX) { G.p.KdMAX++; }
            mutated = true;
            driver_mutated = true;
        }
        // passenger mutation
        if((G.rn.Double() <(G.p.Tp * G.p.mu))) {
            kp++;
            if (kp > G.p.KpMAX) { G.p.KpMAX++; }
            mutated = true;
        }

        if (mutated) {
            parentID = progenyID;
            progenyID = G.p.progenyNextID;
            G.progenyToParentIDs[progenyID] = parentID;
            G.driver_status[progenyID] = (driver_mutated) ? kd : 0; // ryan
            G.p.progenyNextID++;
        }

        return this;
    }
    Cell3D Divide(){
        //int nDivOptions=G.HoodToIs(G.neighborhood,G.divIs,Xsq(),Ysq(),Zsq());
        int nDivOptions = G.MapEmptyHood(G.neighborhood,Xsq(),Ysq(),Zsq());

        int nextAgentID;
        int actualNDivOptions = 0;

        for (int i = 0; i < nDivOptions; i++) {

            //nextAgentID = G.divIs[i];
            nextAgentID = G.neighborhood[i];

            // if invasion hasn't occurred, remain constrained to the ducts
            if (!G.p.INVASION) {
                if ((G.GetAgent(nextAgentID)==null) && (G.ductGrid.Get(nextAgentID) == G.p.OPEN_DUCT)) {
                    //G.divIs[actualNDivOptions] = nextAgentID;
                    G.neighborhood[actualNDivOptions] = nextAgentID;
                    actualNDivOptions ++;
                }
            } else {
                if ((G.GetAgent(nextAgentID)==null)) {
                    //G.divIs[actualNDivOptions] = nextAgentID;
                    G.neighborhood[actualNDivOptions] = nextAgentID;
                    actualNDivOptions ++;
                }
            }


        }

        if(actualNDivOptions==0){ return null; }

        //nextAgentID = G.divIs[G.rn.Int(actualNDivOptions)];
        nextAgentID = G.neighborhood[G.rn.Int(actualNDivOptions)];

        return G.NewAgentSQ(nextAgentID).Init(this.kp,this.kd,this.progenyID,this.parentID).Mutate();
    }
    Cell3D DivideNoMutate(){

        int nDivOptions = G.MapEmptyHood(G.neighborhood,Xsq(),Ysq(),Zsq());
        int nextAgentID;
        int actualNDivOptions = 0;

        // iterate over nDivOptions neighbors, see if they're inside ductal network
        for (int i = 0; i < nDivOptions; i++) {
            nextAgentID = G.neighborhood[i];
            if ((G.GetAgent(nextAgentID)==null) && (G.ductGrid.Get(nextAgentID) == G.p.OPEN_DUCT)) {
                //G.divIs[actualNDivOptions] = nextAgentID;
                G.neighborhood[actualNDivOptions] = nextAgentID;
                actualNDivOptions ++;
            }
        }

        if(actualNDivOptions==0){ return null; }
        nextAgentID = G.neighborhood[G.rn.Int(actualNDivOptions)];

        return G.NewAgentSQ(nextAgentID).Init(this.kp,this.kd,this.progenyID,this.parentID);
    }
    void Step(){
        if (G.p.mode == 1) {
            // mode # 1: P-D ON BIRTH
            if(G.rn.Double()<(Math.pow(1.0+G.p.sd,(double)kd)/Math.pow(1.0+G.p.sp,(double)kp)*G.p.BIRTH_RATE)){
                Divide();
            }

            // check if death event
            if(G.rn.Double()<(G.p.DEATH_RATE )){
                Dispose();
                return;
            }
        } else {
            // mode # 2: P ON BIRTH ; D ON DEATH
            if(G.rn.Double()<(1.0/Math.pow(1.0+G.p.sp,(double)kp)*G.p.BIRTH_RATE)){
                Divide();
            }

            // check if death event
            if(G.rn.Double()<(G.p.DEATH_RATE / Math.pow(1.0+G.p.sd,(double)kd))){
                Dispose();
                return;
            }
        }
    }
    void StepNoDeath(){
        if (G.p.mode == 1) {
            // mode # 1: P-D ON BIRTH
            if(G.rn.Double()<(Math.pow(1.0+G.p.sd,(double)kd)/Math.pow(1.0+G.p.sp,(double)kp)*G.p.BIRTH_RATE)){
                DivideNoMutate();
            }
        } else {
            // mode # 2: P ON BIRTH ; D ON DEATH
            if(G.rn.Double()<(1.0/Math.pow(1.0+G.p.sp,(double)kp)*G.p.BIRTH_RATE)){
                DivideNoMutate();
            }
        }
    }

}

public class PD3D extends AgentGrid3D<Cell3D> {

    Parameters3D p = new Parameters3D();

    // these aren't in param class
    public double[] finalValues= new double[8];
    public int[] progenyToParentIDs = new int[p.progenyArrayLength]; //7500000
    public int[] driver_status = new int[p.progenyArrayLength]; // ryan


    public ArrayList<int[]> duct;
    Grid3Dint ductGrid = new Grid3Dint(p.sideX,p.sideY,p.sideZ);
    int[]neighborhood = MooreHood3D(false);
    Rand rn=new Rand(1);

    PD3D(int sideX, int sideY, int sideZ, Parameters3D params, String ductFilename){
        super(sideX,sideY, sideZ, Cell3D.class, false, false, false);

        this.p = params;
        this.p.KdMAX = params.kd0; // has to match kd0;

        // read in boundaries (ductal structure)
        FileIO reader = new FileIO(ductFilename,"r"); // no_structure, duct_structure.csv
        this.duct = reader.ReadInts(",");
        this.buildDuctGrid();

        // initializer of cancer cells
        int totalCells = 0;

        // start at prescribed z level
        int z = 10;
        boolean foundCell = false;
        while (!foundCell) {
            int x = rn.Int(sideX);
            int y = rn.Int(sideY);
            if ((ductGrid.Get(x,y,z) == 1)) {
                NewAgentSQ(x,y,z).Init(0,this.p.kd0,1,0);
                totalCells++;
                foundCell = true;
            }
        }

        // this "grows" the tumor to initialCells size (all identical cells / no mutations)
        while (totalCells < this.p.initialCells) {
            OriginalStepNoDeath(this.p.initialCells);
            totalCells = Pop();
        }
    }

    // constructor used to draw ducts only
    PD3D(Parameters3D parameters3D, String ductFilename){
        super(parameters3D.sideX, parameters3D.sideY, parameters3D.sideZ, Cell3D.class, false, false, false);

        // read in boundaries (skin or ductal structure)
        FileIO reader = new FileIO(ductFilename,"r"); // no_structure, duct_structure.csv
        this.duct = reader.ReadInts(",");
        this.buildDuctGrid();
    }


    void buildDuctGrid() {
        for (int x = 0; x < ductGrid.xDim; x++) {
            for (int y = 0; y < ductGrid.yDim; y++) {
                for (int z = 0; z < ductGrid.zDim; z++) {
                    this.ductGrid.Set(x,y,z,this.duct.get(z*ductGrid.xDim + x)[y]);
                }
            }
        }
    }

    void OriginalStepNoDeath(int initialCells){
        for (Cell3D c:this) {
            c.StepNoDeath();
            if (this.Pop() >= initialCells) { break; }
        }
        CleanShuffle(rn);
    }

    void OriginalStep(){
        for (Cell3D c:this) {
            c.Step();
        }
        CleanShuffle(rn);
    }


    /*
        SingleSim()

            - runs a single simulation for parameters3D.totalTime time steps
            - invasion occurs after <kd> = 2
                - this happens by setting bigModel.p.INVASION = true;
            - outputs Muller information, basic population statistics ("everything" arrays) and GIFs
     */

    public static void SingleSim(Parameters3D parameters3D, String foldername, String ductFilename) {


        int sideX = parameters3D.sideX;
        int sideY = parameters3D.sideY;
        int sideZ = parameters3D.sideZ;

        // single big model
        PD3D bigModel = new PD3D(sideX, sideY, sideZ, parameters3D, ductFilename);
        int expectedProgeny = bigModel.progenyToParentIDs.length;


        bigModel.p.INVASION = false;

        System.out.println(bigModel.Pop());

        // OUTPUT FILES
        String baseFilename =  foldername + "everything_" + "br" + Integer.toString((int) (parameters3D.BIRTH_RATE*100)) + "dr" + Integer.toString((int) (parameters3D.DEATH_RATE*100)) + "sp" + Integer.toString((int) (parameters3D.sp * 100000)) + "_sd" + Integer.toString((int) (parameters3D.sd * 100000)) + "_Tp" + Integer.toString((int) (parameters3D.Tp )) + "_Td" + Integer.toString((int) (parameters3D.Td )) + "_mode" + Integer.toString(parameters3D.mode);
        FileIO everythingBigOutputFile = new FileIO((baseFilename + "_" + "Big.csv"), "w");
        double[][] everythingBig = new double[8][parameters3D.totalTime / parameters3D.modifier];
        int[][] mullerDriversBig = new int[20][parameters3D.totalTime / parameters3D.modifier];   // track the first 20 drivers (no divisions)
        int[][] mullerGeneticBig = new int[expectedProgeny][parameters3D.totalTime / parameters3D.modifier];

        // initialize all arrays
        for (int jjj = 0; jjj < (parameters3D.totalTime / parameters3D.modifier); jjj++) {
            for (int iii = 0; iii < 8; iii++) {
                everythingBig[iii][jjj] = (iii == 0) ? (double) parameters3D.modifier * jjj : 0.0;
            }
            for (int iii = 0; iii < 20; iii++) {
                mullerDriversBig[iii][jjj] = 0;
            }
            for (int iii = 0; iii < expectedProgeny; iii++) {
                mullerGeneticBig[iii][jjj] = 0;
            }
        }

        OpenGL3DWindow vis=new OpenGL3DWindow("TumorVis", 1000,1000,sideX,sideY,sideZ);

        // draw Z
        String ZgifFilename = foldername + "Z.gif";
        String XgifFilename = foldername + "X.gif";
        String YgifFilename = foldername + "Y.gif";

        GifMaker ZGif = new GifMaker(ZgifFilename, 100,true);
        GifMaker XGif = new GifMaker(XgifFilename, 100,true);
        GifMaker YGif = new GifMaker(YgifFilename, 100,true);

        // wins
        UIWindow Zwin = new UIWindow("Z Gif", true);
        UIWindow Xwin = new UIWindow("X Gif", true);
        UIWindow Ywin = new UIWindow("Y Gif", true);

        // grids
        UIGrid ZVis = new UIGrid(sideX, sideY, 1);
        UIGrid XVis = new UIGrid(sideY, sideZ, 1);
        UIGrid YVis = new UIGrid(sideX, sideZ, 1);

        //guis
        Xwin.AddCol(0, XVis);
        Xwin.RunGui();
        Ywin.AddCol(0, YVis);
        Ywin.RunGui();
        Zwin.AddCol(0, ZVis);
        Zwin.RunGui();



        System.out.println("set up GUIs");

        int i = 0;
        int j = 0;
        int modifier_gif = parameters3D.modifier;

        while (true) {

            if (i < parameters3D.totalTime) {

                // check average driver, then invade
                // pre-invasion: tumor cells constrained inside duct
                // post-invasion: tumor cells have no constraints for which grid points they occupy
                if (AverageDriverNumber(bigModel) > 2.0) {
                    bigModel.p.INVASION = true;
                }

                // make some GIFs
                if (i % modifier_gif==0) {

                    // draw tumor cells
                    // unfortunately, the OpenGL3DWindow has no GIF capabilities, so just print out PNGs:
                    SurfaceDraw(vis, sideX,sideY,sideZ, bigModel);
                    vis.ToPNG(foldername + "img" + Integer.toString(i) + ".png");

                    DrawZ(ZVis,sideX,sideY,sideZ,bigModel);
                    ZGif.AddFrame(ZVis);

                    DrawX(XVis,sideX,sideY,sideZ,bigModel);
                    XGif.AddFrame(XVis);

                    DrawY(YVis,sideX,sideY,sideZ,bigModel);
                    YGif.AddFrame(YVis);


                    // write out
                    Output3DMatrix(bigModel, (i / modifier_gif), foldername);

                }

                // WRITE OUT VALUES
                if (i % parameters3D.modifier == 0) {

                    if (bigModel.Pop() < 10) { break; }
                    System.out.println("time: " + i);

                    // SAVE everythingBig
                    everythingBig[1][j] = bigModel.Pop();               // total pop size
                    everythingBig[2][j] = GetDiversity(bigModel, 1);    // pass div
                    everythingBig[3][j] = GetDiversity(bigModel, 0);    // driv div
                    everythingBig[4][j] = CountMutators(bigModel, 1);   // pass num
                    everythingBig[5][j] = CountMutators(bigModel, 0);   // drivers num
                    everythingBig[6][j] = bigModel.p.KpMAX;             // max pass
                    everythingBig[7][j] = bigModel.p.KdMAX;             // max driv

                    for (int k = 0; k < (sideX * sideY * sideZ); k++) {
                        Cell3D c = bigModel.GetAgent(k);
                        if (c != null) {
                            mullerDriversBig[c.kd][j] = mullerDriversBig[c.kd][j] + 1;
                            mullerGeneticBig[c.progenyID][j] = mullerGeneticBig[c.progenyID][j] + 1;
                        }
                    }

                    // this is kind of a hack to make sure the array is big enough:
                    System.out.println("Next progeny ID is: " + bigModel.p.progenyNextID);
                    System.out.println(bigModel.Pop());
                    j++;
                }

                bigModel.OriginalStep();
                i++;
            } else {

                // OUTPUT "EVERYTHING" ARRAYS TO FILE
                BuildEverything(everythingBigOutputFile, everythingBig, parameters3D);

                // OUTPUT DRIVER MULLER PLOTS:
                FileIO mullerBigOutputFile = new FileIO((foldername + "driverClonesBig.csv"), "w");
                WriteDriverClones(mullerBigOutputFile,mullerDriversBig, parameters3D);

                System.out.println("before reduce");

                ////// genetic clones of REDUCED, BIG
                FileIO bigParentsReduced = new FileIO((foldername + "parentsBig.csv"), "w");
                FileIO bigGeneticMullerReduced = new FileIO((foldername + "clonesBig.csv"), "w");
                FileIO bigDriverStatusReduced = new FileIO((foldername + "driverStatusBig.csv"), "w");
                ReduceParents(bigParentsReduced, bigGeneticMullerReduced, bigDriverStatusReduced, mullerGeneticBig, bigModel.progenyToParentIDs, bigModel.driver_status, expectedProgeny, parameters3D);

                System.out.println("after reduce");

                ZGif.Close();
                XGif.Close();
                YGif.Close();
                //vis.Close();

                System.out.println("okay to close");

                break;
            }
        }

        return;
    }

    public static void main(String[] args) {

        // run a single simulation of 3d tumor growth constrained to ductal network
        // tumor invasion occurs when average num. drivers is <kd> = 2

        String folderName = "data-output/";
        String ductFilename = "Model/slices/duct_structure3.csv";
        Parameters3D parameters3D = new Parameters3D();
        parameters3D.totalTime = 300; // was 1000
        parameters3D.modifier = 20;
        parameters3D.initialCells = 500;

        SingleSim(parameters3D, folderName, ductFilename);

        return;
    }

    /*

        DrawX()
            - two dimensional side view (a projection, not a slice) from x-camera POV

     */

    public static void DrawX(UIGrid vis, int sideX, int sideY, int sideZ, PD3D model){
        for (int yy = 0; yy < sideY; yy++) {
            for (int zz = 0; zz < sideZ; zz++) {
                for (int xx = 0; xx < sideX; xx++) {
                    Cell3D c=model.GetAgent(xx,yy,zz);

                    if(c!=null){
                        int color = CategorialColor((c.kd - 1) % 19);// + rn.nextInt(1000) - 500;
                        vis.SetPix(c.Ysq(), c.Zsq(),color);
                        break;
                    } else {
                        vis.SetPix(yy, zz,RGB(1,1,1));
                    }
                }
            }
        }
    }

    /*

        DrawY()
            - two dimensional side view (a projection, not a slice) from y-camera POV

     */

    public static void DrawY(UIGrid vis, int sideX, int sideY, int sideZ, PD3D model){
        for (int xx = 0; xx < sideX; xx++) {
            for (int zz = 0; zz < sideZ; zz++) {
                for (int yy = 0; yy < sideY; yy++) {
                    Cell3D c=model.GetAgent(xx,yy,zz);
                    if(c!=null){
                        int color = CategorialColor((c.kd - 1) % 19);// + rn.nextInt(1000) - 500;
                        vis.SetPix(c.Xsq(), c.Zsq(),color);
                        break;
                    } else {
                        vis.SetPix(xx, zz,RGB(1,1,1));
                    }
                }
            }
        }
    }


    /*

        DrawZ()
            - two dimensional side view (a projection, not a slice) from z-camera POV

     */

    public static void DrawZ(UIGrid vis, int sideX, int sideY, int sideZ, PD3D model){

        System.out.println("pop is : " + model.Pop());
        int testing = 0;
        // from z = 0 side
        for (int xx = 0; xx < sideX; xx++) {
            for (int yy = 0; yy < sideY; yy++) {
                for (int zz = 0; zz < sideZ; zz++) {
                    Cell3D c=model.GetAgent(xx,yy,zz);
                    if(c!=null){
                        testing++;
                        int color = CategorialColor((c.kd - 1) % 19);// + rn.nextInt(1000) - 500;
                        vis.SetPix(c.Xsq(), c.Ysq(),color);
                        break;
                    } else {
                        vis.SetPix(xx, yy,RGB(1,1,1));
                    }
                }
            }
        }
        System.out.println("testing is " + testing  );
    }


    /*

        SurfaceDraw()
            - Uses OpenGL3DWindow to display a three dimensional rendering of all tumor cells
            - shaded from dark (back, left, bottom corner) to light (top, right, front corner)

     */

    public static void SurfaceDraw(OpenGL3DWindow vis, int sideX, int sideY, int sideZ, PD3D model){
        vis.ClearBox(RGB(1,1,1),RGB(0,0,0));//used to clear gui

        int maxZ = 0;
        for (int zz = 0; zz < sideZ; zz++) {
            // check this z layer, if cell
            double Zdepth = (double) zz / sideZ; // draw from back (use 1  - to draw from front)

            for (int xx = 0; xx < sideX; xx++) {
                for (int yy = 0; yy < sideY; yy++) {

                    Cell3D c=model.GetAgent(xx,yy,zz);
                    double YDepth = (double) yy / sideY;
                    double totalDepth = (Zdepth + YDepth) / 2.0;
                    double frac = 0.9;

                    if(c!=null){
                        int color = RGB(totalDepth*frac,totalDepth*frac,totalDepth*frac);//RGB(depth,0,depth2);
                        vis.Circle(c.Xsq(),c.Ysq(),sideZ - c.Zsq(),0.7, color); // this draws from the back
                        maxZ = zz;
                    }
                }
            }
        }

        vis.Update();
    }

    /*

        DrawDuctOnly()
            - Renders a three dimensional drawing of ductal network structure only (no tumor)

     */

    public static void DrawDuctOnly(OpenGL3DWindow vis, Parameters3D parameters3D, PD3D model){
        int maxZ = parameters3D.sideZ;
        double circleSize = 2;
        int modifier = 2;

        for (int zz = 0; zz < maxZ; zz+=modifier) {
            double Zdepth = (double) zz / parameters3D.sideZ;
            for (int xx = 0; xx < parameters3D.sideX; xx+=modifier) {
                for (int yy = 0; yy < parameters3D.sideY; yy+=modifier) {
                    if (model.ductGrid.Get(xx, yy, zz) == 1) {
                        vis.Circle(xx, yy, parameters3D.sideZ - zz, circleSize, RGB(Zdepth, 0, 1 - Zdepth));
                    }
                }
            }
        }

        vis.Update();
    }

    /*

        GetDiversity()
            - iterates through all tumor cells to find Shannon entropy, H
            - H = exp(- sum( p_i  * log (p_i) ) )
            - where p_i = fraction of cells in tumor population with i driver (or passengers)
            - use passBool = 0 to get driver diversity
            - use passBool = 1 to get passenger diversity

     */

    public static double GetDiversity(PD3D model, int passBool) {
        if (passBool == 0) {
            // count drivers!
            int kVec[] = new int[model.p.KdMAX + 1];
            for (Cell3D c : model) {
                kVec[c.kd - 1]++;
            }

            // calculate diversity
            double sum = 0;
            for (int i = 0; i <= model.p.KdMAX; i++) {
                if (kVec[i] > 0) {
                    sum += (double)((double)kVec[i] / (double)model.Pop())*Math.log10((double)kVec[i] / (double)model.Pop());
                }
            }
            return Math.exp(- sum);

        } else {
            // count passengers!
            int kVec[] = new int[model.p.KpMAX + 1];
            for (Cell3D c : model) {
                kVec[c.kp]++;
            }

            // calculate diversity
            double sum = 0;
            for (int i = 0; i <= model.p.KpMAX; i++) {
                if (kVec[i] > 0) {
                    sum += (double)((double)kVec[i] / (double)model.Pop())*Math.log10((double)kVec[i] / (double)model.Pop());
                }
            }
            return Math.exp(- sum);

        }
    }

    /*

        AverageDriverNumber()
            - iterates through all tumor cells to find average driver number, per cell: <k>
            - simulation begins with all tumor cells k_d = 1 so <k> = 1
            - <k> = sum (k_d) / N where N is tumor population size

     */

    public static double AverageDriverNumber(PD3D model) {
        int totalDrivers = CountMutators(model, 0);
        return (double)totalDrivers / (double)model.Pop();
    }

    /*

        CountMutators()
            - iterates through all tumor cells to sum k_d or k_p for all cells
            - use passBool = 0 to get driver mutations
            - use passBool = 1 to get passenger mutations

     */

    public static int CountMutators(PD3D model, int passBool) {
        int sum = 0;
        if (passBool == 0) {
            // count drivers!

            for (Cell3D c : model) {
                sum += c.kd;
            }
            return sum;

        } else {
            // count passengers!
            for (Cell3D c : model) {
                sum += c.kp;
            }

            return sum;
        }


    }


    /*

        BuildEverything()
            - outputs a file of all important tumor metrics
            - rows are metrics; columns are timepoints
            - Metrics:
                - Row 0: time
                - Row 1: population size, N
                - Row 2: Shannon diversity of passengers, H_p
                - Row 3: Shannon diversity of drivers, H_d
                - Row 4: Total number of passenger mutations - see CountMutators()
                - Row 5: Total number of driver mutations - see CountMutators()
                - Row 6: Maximum passenger mutation found in ANY cell in population
                - Row 7: Maximum driver mutation found in ANY cell in population

     */

    public static void BuildEverything(FileIO everythingBigOutputFile, double[][] everythingBig, Parameters3D parameters3D) {
        StringBuilder sb = new StringBuilder();

        for (int jj = 0; jj < 8; jj++) {
            sb = new StringBuilder();
            for (int ii = 0; ii < (parameters3D.totalTime / parameters3D.modifier) - 1; ii++) { sb.append(everythingBig[jj][ii] + ","); }
            sb.append(everythingBig[jj][(parameters3D.totalTime / parameters3D.modifier) - 1] + "\n");
            everythingBigOutputFile.Write(sb.toString());
        }

        everythingBigOutputFile.Close();
    }

    /*

        WriteDriverClones()
            - outputs a file where rows are individual (unique) clone populations

     */

    public static void WriteDriverClones(FileIO mullerBigOutputFile, int[][] mullerDriversBig, Parameters3D parameters3D) {
        // first, output time:
        StringBuilder sb = new StringBuilder();

        for (int time = 0; time < (parameters3D.totalTime / parameters3D.modifier) - 1; time++) { sb.append(time* parameters3D.modifier + ","); }
        sb.append(parameters3D.modifier*((parameters3D.totalTime / parameters3D.modifier) - 1) + "\n");
        mullerBigOutputFile.Write(sb.toString());

        for (int jj = 0; jj < 20; jj++) {
            sb = new StringBuilder();
            for (int ii = 0; ii < (parameters3D.totalTime / parameters3D.modifier) - 1; ii++) {
                sb.append(mullerDriversBig[jj][ii] + ",");
            }
            sb.append(mullerDriversBig[jj][(parameters3D.totalTime / parameters3D.modifier) - 1] + "\n");
            mullerBigOutputFile.Write(sb.toString());
        }

        mullerBigOutputFile.Close();

    }

    /*

        ReduceParents()
            - removes clones that are dead and gone
            - reconnects parents together so that there is a continuous lineage when parent dies and child lives on
            - writes everything out

     */


    public static void ReduceParents(FileIO parentsReduced, FileIO geneticMullerReduced, FileIO driverStatusReduced, int[][] mullerGenetic, int[] progenyToParentIDs, int[] driver_status, int expectedProgeny, Parameters3D parameters3D) {
        int rowmax = 0;
        ArrayList<Integer> reducedParents = new ArrayList<>();
        ArrayList<Integer> reducedProgeny = new ArrayList<>();
        ArrayList<Integer> reducedDriverStatus = new ArrayList<>();
        reducedParents.add(new Integer(-1)); // add 0th parent
        reducedProgeny.add(new Integer(0)); // add 0th progeny


        for (int jj = 1; jj < expectedProgeny; jj++) {
            rowmax = 0;
            for (int ii = 0; ii < (parameters3D.totalTime / parameters3D.modifier) - 1; ii++) {
                rowmax = (mullerGenetic[jj][ii] > rowmax) ? mullerGenetic[jj][ii]: rowmax;
            }

            if (rowmax > parameters3D.delete_thresh) {

                int mySupposedParent = progenyToParentIDs[jj];

                // check if my parent is in there:
                while (!reducedProgeny.contains(new Integer(mySupposedParent))) {

                    // parent is parent of my parent, lol -- is that parent in there?
                    mySupposedParent = progenyToParentIDs[progenyToParentIDs[mySupposedParent]];

                    // breaks if it is in there, eventually up the chain
                }

                // if my parent is in the reduced thing, add me.
                //reducedParentsBig.add(new Integer(mySupposedParent));
                reducedParents.add(new Integer(reducedProgeny.indexOf(mySupposedParent)));
                reducedProgeny.add(new Integer(jj));
                reducedDriverStatus.add(new Integer(driver_status[jj])); // ryan

            }
        }


        // first, output time:
        StringBuilder sb = new StringBuilder();
        StringBuilder sb2 = new StringBuilder();
        StringBuilder sb3 = new StringBuilder(); //ryan
        for (int time = 0; time < (parameters3D.totalTime / parameters3D.modifier) - 1; time++) { sb.append(time* parameters3D.modifier + ","); }
        sb.append(parameters3D.modifier*((parameters3D.totalTime / parameters3D.modifier) - 1) + "\n");
        geneticMullerReduced.Write(sb.toString()); // file name


        // out gene muller big to file
        int ii_reduced = 0;
        for (int jj = 1; jj < expectedProgeny; jj++) {
            sb = new StringBuilder();

            // check if it's in reduced
            if (reducedProgeny.contains(new Integer(jj))) {
                for (int ii = 0; ii < (parameters3D.totalTime / parameters3D.modifier) - 1; ii++) {
                    sb.append(mullerGenetic[jj][ii] + ",");
                }
                sb.append(mullerGenetic[jj][(parameters3D.totalTime / parameters3D.modifier) - 1] + "\n");

                // output string, including newline \n
                geneticMullerReduced.Write(sb.toString());

                // also output this particular parent

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




    public static void Output3DMatrix(PD3D model, int j, String foldername) {
        StringBuilder sb = new StringBuilder();
        String baseFilename =  foldername + "matrix" + Integer.toString(j);
        FileIO outputFile = new FileIO((baseFilename + ".csv"), "w");

        System.out.println("xDim : " + model.xDim);
        System.out.println("yDim : " + model.yDim);
        System.out.println("zDim : " + model.zDim);

        for (int zz = 0; zz < model.zDim; zz++) {
            for (int xx = 0; xx < model.xDim; xx++) {
                for (int yy = 0; yy < (model.yDim-1); yy++) {
                    Cell3D c = model.GetAgent(xx, yy, zz);
                    sb.append((c != null) ? 1 + "," : 0 + ",");

                }

                // last y dim & new line
                Cell3D c = model.GetAgent(xx, model.yDim-1, zz);
                sb.append((c != null) ? 1 + "\n" : 0 + "\n");

            }
        }

        outputFile.Write(sb.toString());

    }


}
