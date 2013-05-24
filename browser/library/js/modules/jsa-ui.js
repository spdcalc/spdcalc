define(
    [
        'jquery',
        'stapes',
        'phasematch',
        'modules/heat-map',
        'modules/line-plot',
        'tpl!templates/jsa-layout.tpl'
    ],
    function(
        $,
        Stapes,
        PhaseMatch,
        HeatMap,
        LinePlot,
        tplJSALayout
    ) {

        'use strict';


        var con = PhaseMatch.constants;
        var defaults = {
            
        };

        /**
         * @module JSAUI
         * @implements {Stapes}
         */
        var jsaUI = Stapes.subclass({

            /**
             * Mediator Constructor
             * @return {void}
             */
            constructor : function( config ){

                var self = this;

                self.options = $.extend({}, defaults, config);

                self.initPhysics();

                self.el = $( tplJSALayout.render() );

                var margins = {
                    top: 60,
                    right: 40,
                    left: 80,
                    bottom: 60
                };

                // init plot
                self.plot = new HeatMap({
                    el: self.el.find('.heat-map-wrapper').get( 0 ),
                    margins: margins,
                    width: 480,
                    height: 480,
                    labels: {
                        x: 'Wavelength of Signal (nm)',
                        y: 'Wavelength of Idler (nm)'
                    },
                    xrange: [ 0, 200 ],
                    yrange: [ 0, 100 ],
                    format: '.0f'
                });
                self.plot.setTitle('Joint spectral amplitude');
                self.elPlot = $(self.plot.el);

            },

            initPhysics: function(){

                // initialize physics if needed...
            },

            /**
             * Connect to main app
             * @return {void}
             */
            connect : function( app ){

                var self = this
                    ;

                self.parameters = app.parameters;

                // connect to the app events
                app.on({

                    calculate: self.refresh

                }, self);

                // auto draw
                self.refresh();
                
            },

            disconnect: function( app ){

                // disconnect from app events
                app.off( 'calculate', self.refresh );
            },

            resize: function(){

                var self = this
                    ,par = self.elPlot.parent()
                    ,width = par.width()
                    ,height = $(window).height()
                    ,dim = Math.min( width, height ) - 100 // - margin
                    ;

                // if (dim > 400){ 
                //     dim = 400;
                // }

                // self.plot.resize( dim, dim );
                self.draw();
            },

            getMainPanel: function(){
                return this.el;
            },

            refresh: function(){

                var self = this;
                self.calc( self.parameters.getProps() );
                self.draw();
            },

            calc: function( props ){

                // @TODO: move this to a control bar
                props.lambda_i = 1/(1/props.lambda_p - 1/props.lambda_s);
                var dim = 200;
                // var l_start = 1500 * con.nm;
                // var l_stop = 1600 * con.nm; 
                var threshold = 0.5;
                var lsi = PhaseMatch.autorange_lambda(props, threshold);
                var l_start = Math.min(lsi[0], lsi[1]);
                var l_stop =  Math.max(lsi[0], lsi[1]);
                // console.log("max, min ",threshold,  l_start/1e-9, l_stop/1e-9);
                var data1d = [];

                var self = this
                    ,PM = PhaseMatch.calc_JSA(
                        props, 
                        l_start, 
                        l_stop, 
                        lsi[2],
                        lsi[3], 
                        dim
                    )
                    ;

                self.data = PM;
                // multiplication is faster than division
                self.plot.setXRange([l_start * 1e9, l_stop * 1e9]);
                self.plot.setYRange([lsi[2] * 1e9, lsi[3] * 1e9]);

                

            },

            draw: function(){

                var self = this
                    ,data = self.data
                    ;

                if (!data){
                    return this;
                }

                self.plot.plotData( data );

               
            }
        });

        return function( config ){

            return new jsaUI( config );
        };
    }
);