/**
 * These are the properties that are used to calculate phasematching
 */

(function(){

    // /**
    //  * Rotation object
    //  */
    // var Rotation = function(){

    //     this.Sx = 0;
    //     this.Sy = 0;
    //     this.Sz = 0;
    // };

    // Rotation.prototype = {

    //     set: function( theta, phi, theta_s, phi_s ){

    //         // First get the ransfomration to lambda_p coordinates
    //         var S_x = Math.sin(theta_s)*Math.cos(phi_s);
    //         var S_y = Math.sin(theta_s)*Math.sin(phi_s);
    //         var S_z = Math.cos(theta_s);

    //         // Transform from the lambda_p coordinates to crystal coordinates
    //         var SR_x = Math.cos(theta)*Math.cos(phi)*S_x - Math.sin(phi)*S_y + Math.sin(theta)*Math.cos(phi)*S_z;
    //         var SR_y = Math.cos(theta)*Math.sin(phi)*S_x + Math.cos(phi)*S_y + Math.sin(theta)*Math.sin(phi)*S_z;
    //         var SR_z = -Math.sin(theta)*S_x  + Math.cos(theta)*S_z;
            
    //         // Normalambda_ize the unit vector
    //         // FIX ME: When theta = 0, Norm goes to infinity. This messes up the rest of the calculations. In this
    //         // case I think the correct behaviour is for Norm = 1 ?
    //         var Norm =  Math.sqrt(sq(S_x) + sq(S_y) + sq(S_z));
    //         this.Sx = SR_x/(Norm);
    //         this.Sy = SR_y/(Norm);
    //         this.Sz = SR_z/(Norm);
    //     }
    // };

    // PhaseMatch.Rotation = Rotation;

    /**
     * SPDCprop
     */
    var SPDCprop = function(){
        this.init();

    };

    SPDCprop.prototype = {
        init:function(){
            var con = PhaseMatch.constants;
            this.lambda_p = 775 * con.nm;
            this.lambda_s = 1550 * con.nm;
            this.lambda_i = 1550 * con.nm;
            this.Types = ["o -> o + o", "e -> o + o", "e -> e + o", "e -> o + e"];
            this.Type = this.Types[1];
            this.theta = 19.8371104525 *Math.PI / 180;
            // this.theta = 19.2371104525 *Math.PI / 180;
            this.phi = 0;
            this.theta_s = 0; // * Math.PI / 180;
            this.theta_i = 0;
            this.phi_s = 0;
            this.phi_i = 0;
            this.poling_period = 1000000;
            this.L = 20000 * con.um;
            this.W = 500 * con.um;
            this.p_bw = 1;
            this.phase = false;
            this.apodization = 1;
            this.apodization_FWHM = 1000 * con.um;
            this.crystal = new PhaseMatch.BBO();
            //Other functions that do not need to be included in the default init
            this.S_p = this.calc_Coordinate_Transform(this.theta, this.phi, 0, 0);
            this.S_s = this.calc_Coordinate_Transform(this.theta, this.phi, this.theta_s, this.phi_s);
            this.S_i = this.calc_Coordinate_Transform(this.theta, this.phi, this.theta_i, this.phi_i);

            this.n_p = this.calc_Index_PMType(this.lambda_p, this.Type, this.S_p, "pump");
            this.n_s = this.calc_Index_PMType(this.lambda_s, this.Type, this.S_s, "signal");
            this.n_i = this.calc_Index_PMType(this.lambda_i, this.Type, this.S_i, "idler");
        },
            // this.autocalcTheta = false;
            // this.calc_theta= function(){
            //     //unconstrained minimization
            //     if this.autocalcTheta{}
            //     return this.theta = answer
            // }
        calc_Coordinate_Transform : function (theta, phi, theta_s, phi_s){
            //Should save some calculation time by defining these variables.
            var SIN_THETA = Math.sin(theta);
            var COS_THETA = Math.cos(theta);
            var SIN_THETA_S = Math.sin(theta_s);
            var COS_THETA_S = Math.cos(theta_s);
            var SIN_PHI = Math.sin(phi);
            var COS_PHI = Math.cos(phi);
            var SIN_PHI_S = Math.sin(phi_s);
            var COS_PHI_S = Math.cos(phi_s);


            var S_x = SIN_THETA_S*COS_PHI_S;
            var S_y = SIN_THETA_S*SIN_PHI_S;
            var S_z = COS_THETA_S;

            // Transform from the lambda_p coordinates to crystal coordinates
            var SR_x = COS_THETA*COS_PHI*S_x - SIN_PHI*S_y + SIN_THETA*COS_PHI*S_z;
            var SR_y = COS_THETA*SIN_PHI*S_x + COS_PHI*S_y + SIN_THETA*SIN_PHI*S_z;
            var SR_z = -SIN_THETA*S_x  + COS_THETA*S_z;
            
            // Normalambda_ize the unit vector
            // @TODO: When theta = 0, Norm goes to infinity. This messes up the rest of the calculations. In this
            // case I think the correct behaviour is for Norm = 1 ?
            var Norm =  Math.sqrt(sq(S_x) + sq(S_y) + sq(S_z));
            var Sx = SR_x/(Norm);
            var Sy = SR_y/(Norm);
            var Sz = SR_z/(Norm);

            return [Sx, Sy, Sz];
        },

        calc_Index_PMType : function (lambda, Type, S, photon){
            var ind = this.crystal.indicies(lambda);

            var nx = ind[0];
            var ny = ind[1];
            var nz = ind[2];

            var Sx = S[0];
            var Sy = S[1];
            var Sz = S[2];

            var B = sq(Sx) * (1/sq(ny) + 1/sq(nz)) + sq(Sy) *(1/sq(nx) + 1/sq(nz)) + sq(Sz) *(1/sq(nx) + 1/sq(ny));
            var C = sq(Sx) / (sq(ny) * sq(nz)) + sq(Sy) /(sq(nx) * sq(nz)) + sq(Sz) / (sq(nx) * sq(ny));
            var D = sq(B) - 4 * C;

            var nslow = Math.sqrt(2/ (B + Math.sqrt(D)));
            var nfast = Math.sqrt(2/ (B - Math.sqrt(D)));
            //nfast = o, nslow = e

            var n = 1;

            switch (Type){

                case "e -> o + o":
                    if (photon === "pump") { n = nslow;}
                    else { n = nfast;}
                break;
                case "e -> e + o":
                    if (photon === "idler") { n = nfast;}
                    else {n = nslow;}
                break;
                case "e -> o + e":
                    if (photon === "signal") { n = nfast;}
                    else {n = nslow;}
                break;
                default:
                    throw "Error: bad PMType specified";
            }

            return n ;
        },


        set: function( name, val ){


            switch ( name ){

                case 'lambda_p':

                    // this.updateSomethng();
                break;
            }

            this[ name ] = val;
        }
    };

    PhaseMatch.SPDCprop = SPDCprop;
})();

