## =============================================================================
## Cervical Cancer QSP Shiny Dashboard
## Interactive Simulation: Cisplatin CCRT, Pembrolizumab, Bevacizumab, Tisotumab Vedotin
## =============================================================================

library(shiny)
library(shinydashboard)
library(mrgsolve)
library(dplyr)
library(ggplot2)
library(tidyr)
library(DT)
library(plotly)

## ── mrgsolve model code (embedded) ─────────────────────────────────────────
cc_code <- '
$PARAM
CL_CIS=30.0, V1_CIS=15.0, Q_CIS=20.0, V2_CIS=30.0,
CL_PAC=13.2, V1_PAC=6.5, Q_PAC=7.0, V2_PAC=113.0,
CL_BEV=0.207, V1_BEV=2.91, Q_BEV=0.469, V2_BEV=1.91,
CL_PEM=0.213, V1_PEM=3.34, Q_PEM=0.638, V2_PEM=2.68,
CL_TV=0.51, V1_TV=3.35, Q_TV=0.55, V2_TV=2.75, k_dec_MMAE=5.0,
VEGF0=0.20, ksyn_VEGF=1.8, kdeg_VEGF=8.0, kbind_BEV=45.0,
k_adduct=0.02, k_repair=0.10,
alpha_rt=0.30, beta_rt=0.030, k_radiosens=1.6, k_reoxy=0.05, k_repop=0.008,
TV0=40.0, kg=0.010, TV_max=1500.0,
k_kill_RT=1.0, k_kill_Pt=0.005, k_kill_pac=0.003, k_kill_ICI=0.004, k_kill_ADC=0.006, k_kill_bev=0.0015,
SCCAg_0=8.0, ksyn_SCC=0.10, kdeg_SCC=0.25,
HPVload_0=5.0, k_HPV_prod=0.08, k_HPV_clear=0.05, k_HPV_ICI=2.0,
CD8T_0=1.0, k_CD8_in=0.12, k_CD8_out=0.10, k_exhaust=0.35, k_ICI=2.5,
PDL1_0=1.0, k_PDL1_up=0.02, PDL1_max=3.0,
ICI_flag=0, RT_flag=0, CPS_high=1

$CMT CIS_C1 CIS_C2 PAC_C1 PAC_C2 BEV_C1 BEV_C2 PEMBRO_C1 PEMBRO_C2
     TV_ADC_C1 TV_ADC_C2 MMAE_free VEGF Pt_DNA RT_SF TV SCCAg HPVload CD8T PDL1_exp

$MAIN
double eff_ICI = ICI_flag * CPS_high;

$ODE
dxdt_CIS_C1 = -(CL_CIS/V1_CIS)*CIS_C1 - (Q_CIS/V1_CIS)*CIS_C1 + (Q_CIS/V2_CIS)*CIS_C2;
dxdt_CIS_C2 = (Q_CIS/V1_CIS)*CIS_C1 - (Q_CIS/V2_CIS)*CIS_C2;
dxdt_PAC_C1 = -(CL_PAC/V1_PAC)*PAC_C1 - (Q_PAC/V1_PAC)*PAC_C1 + (Q_PAC/V2_PAC)*PAC_C2;
dxdt_PAC_C2 = (Q_PAC/V1_PAC)*PAC_C1 - (Q_PAC/V2_PAC)*PAC_C2;
double BEV_effect = kbind_BEV * BEV_C1 * VEGF;
dxdt_BEV_C1 = -(CL_BEV/V1_BEV)*BEV_C1 - (Q_BEV/V1_BEV)*BEV_C1 + (Q_BEV/V2_BEV)*BEV_C2;
dxdt_BEV_C2 = (Q_BEV/V1_BEV)*BEV_C1 - (Q_BEV/V2_BEV)*BEV_C2;
dxdt_PEMBRO_C1 = -(CL_PEM/V1_PEM)*PEMBRO_C1 - (Q_PEM/V1_PEM)*PEMBRO_C1 + (Q_PEM/V2_PEM)*PEMBRO_C2;
dxdt_PEMBRO_C2 = (Q_PEM/V1_PEM)*PEMBRO_C1 - (Q_PEM/V2_PEM)*PEMBRO_C2;
dxdt_TV_ADC_C1 = -(CL_TV/V1_TV)*TV_ADC_C1 - (Q_TV/V1_TV)*TV_ADC_C1 + (Q_TV/V2_TV)*TV_ADC_C2;
dxdt_TV_ADC_C2 = (Q_TV/V1_TV)*TV_ADC_C1 - (Q_TV/V2_TV)*TV_ADC_C2;
dxdt_MMAE_free = 0.4*(CL_TV/V1_TV)*TV_ADC_C1 - k_dec_MMAE*MMAE_free;
double VEGF_prod = ksyn_VEGF * (TV/TV0);
dxdt_VEGF = VEGF_prod - kdeg_VEGF*VEGF - BEV_effect;
dxdt_Pt_DNA = k_adduct*CIS_C1 - k_repair*Pt_DNA;
double alpha_eff = alpha_rt*(1.0 + (k_radiosens-1.0)*(Pt_DNA/(Pt_DNA+0.3)));
double dose_rate = RT_flag*2.0;
double rt_damage_rate = alpha_eff*dose_rate + beta_rt*dose_rate*dose_rate;
dxdt_RT_SF = rt_damage_rate - k_reoxy*RT_SF;
double ICI_eff = 1.0 + eff_ICI*(k_ICI-1.0);
double exhaustion = k_exhaust*(TV/TV_max)*(PDL1_exp/PDL1_0);
dxdt_CD8T = k_CD8_in*ICI_eff - k_CD8_out*CD8T - exhaustion*CD8T;
dxdt_PDL1_exp = k_PDL1_up*CD8T*(1 - PDL1_exp/PDL1_max);
double grow_term = kg*TV*log(TV_max/TV);
double repop_term = RT_flag*k_repop*TV;
double kill_RT = k_kill_RT*(1-exp(-RT_SF))*TV;
double kill_Pt = k_kill_Pt*Pt_DNA*TV;
double pac_eff = PAC_C1/(PAC_C1+100.0);
double kill_pac = k_kill_pac*pac_eff*TV;
double kill_ICI = k_kill_ICI*eff_ICI*CD8T*TV;
double adc_eff = MMAE_free/(MMAE_free+1.0);
double kill_ADC = k_kill_ADC*adc_eff*TV;
double kill_bev = k_kill_bev*(BEV_C1/(BEV_C1+10.0))*TV;
dxdt_TV = grow_term + repop_term - kill_RT - kill_Pt - kill_pac - kill_ICI - kill_ADC - kill_bev;
double SCC_prod = ksyn_SCC*TV;
dxdt_SCCAg = SCC_prod - kdeg_SCC*SCCAg;
double HPV_clear_eff = k_HPV_clear*(1.0 + eff_ICI*(k_HPV_ICI-1.0));
dxdt_HPVload = k_HPV_prod*(TV/TV0) - HPV_clear_eff*HPVload;

$TABLE
capture CIS_Conc = CIS_C1;
capture PAC_Conc = PAC_C1;
capture BEV_Conc = BEV_C1;
capture PEMBRO_Conc = PEMBRO_C1;
capture TVADC_Conc = TV_ADC_C1;
capture MMAE_lvl = MMAE_free;
capture VEGF_free = VEGF;
capture PtDNA_rel = Pt_DNA;
capture RT_damage = RT_SF;
capture TumorVol = TV;
capture SCCAg_lvl = SCCAg;
capture HPV_rel = HPVload;
capture CD8T_rel = CD8T;
capture PDL1_rel = PDL1_exp;

$INIT
CIS_C1=0,CIS_C2=0,PAC_C1=0,PAC_C2=0,BEV_C1=0,BEV_C2=0,
PEMBRO_C1=0,PEMBRO_C2=0,TV_ADC_C1=0,TV_ADC_C2=0,MMAE_free=0,
VEGF=0.20,TV=40,Pt_DNA=0,RT_SF=0,SCCAg=8.0,HPVload=5.0,CD8T=1.0,PDL1_exp=1.0
'

## Compile once at startup
mod_base <- mcode("cc_shiny", cc_code, quiet=TRUE)

## ── Shiny UI ───────────────────────────────────────────────────────────────
ui <- dashboardPage(
  skin = "red",
  dashboardHeader(title = "Cervical Cancer QSP Dashboard"),

  dashboardSidebar(
    width = 280,
    sidebarMenu(
      menuItem("① Patient profile", tabName="patient",  icon=icon("user")),
      menuItem("② Drug PK",       tabName="pk",       icon=icon("pills")),
      menuItem("③ PD key indicators",   tabName="pd",       icon=icon("flask")),
      menuItem("④ Clinical endpoints", tabName="clinical", icon=icon("chart-line")),
      menuItem("⑤ Scenario comparison", tabName="scenario", icon=icon("layer-group")),
      menuItem("⑥ Biomarkers",    tabName="biomarker", icon=icon("microscope"))
    ),
    hr(),
    h5("── Patient characteristics ──", style="color:#ccc;padding-left:15px"),
    sliderInput("TV0",    "Initial tumour size (cm³)",  5, 200, 40, step=5),
    selectInput("FIGO",   "FIGO 2018 stage",
      choices=c("IB2-IIA (localised)"="early", "IIB-IVA (locally advanced)"="lacc", "IVB/recurrent-metastatic"="rm"), selected="lacc"),
    selectInput("CPS_high","PD-L1 CPS status",
      choices=c("CPS ≥ 1 (positive)"=1, "CPS < 1 (negative)"=0), selected=1),
    sliderInput("SCCAg_0","Initial SCC-Ag (ng/mL)", 1, 50, 8, step=1),
    hr(),
    h5("── Treatment settings ──", style="color:#ccc;padding-left:15px"),
    checkboxInput("add_RT","Concurrent chemoradiotherapy (CCRT)", value=TRUE),
    sliderInput("n_cis","Number of cisplatin doses (weekly)", 3, 6, 6, step=1),
    checkboxInput("add_pembro","Add pembrolizumab (KEYNOTE-A18/826)", value=FALSE),
    checkboxInput("add_bev","Add bevacizumab (GOG-240, R/M)", value=FALSE),
    checkboxInput("add_tv","Tisotumab vedotin (2L+, R/M)", value=FALSE),
    sliderInput("sim_days","Simulation period (days)", 180, 1095, 730, step=30),
    actionButton("run_sim","▶ Run simulation", class="btn-success btn-block")
  ),

  dashboardBody(
    tabItems(
      ## Tab 1: Patient Profile
      tabItem("patient",
        fluidRow(
          valueBoxOutput("box_stage"), valueBoxOutput("box_cps"),
          valueBoxOutput("box_sccag0")
        ),
        fluidRow(
          box(title="Patient characteristics summary", width=6, status="danger",
            tableOutput("patient_summary")
          ),
          box(title="Treatment eligibility matrix", width=6, status="warning",
            tableOutput("eligibility_table")
          )
        ),
        fluidRow(
          box(title="Cervical cancer pathophysiology overview", width=12, status="info",
            HTML('
<table class="table table-bordered table-sm">
<thead><tr style="background:#C0392B;color:white">
  <th>Characteristic</th><th>Locally advanced (LACC)</th><th>Recurrent/metastatic (R/M)</th>
</tr></thead>
<tbody>
<tr><td>Standard treatment</td><td>Cisplatin concurrent chemoradiotherapy (CCRT)</td><td>Platinum-based chemotherapy ± Bev/Pembro</td></tr>
<tr><td>Key pathogenesis</td><td>HPV16/18 E6-p53 / E7-pRb pathway</td><td>Same + selection of treatment-resistant clones</td></tr>
<tr><td>PD-L1 CPS≥1 proportion</td><td>~80-90%</td><td>~80-90%</td></tr>
<tr><td>Key clinical trials</td><td>RTOG-90-01, KEYNOTE-A18</td><td>GOG-240, KEYNOTE-826, innovaTV 301</td></tr>
<tr><td>5-year survival</td><td>IIB: ~65%, IIIB: ~40%</td><td>IVB: ~15-20%</td></tr>
<tr><td>Second-line-plus options</td><td>-</td><td>Tisotumab vedotin (innovaTV 204/301)</td></tr>
</tbody>
</table>'
            )
          )
        )
      ),

      ## Tab 2: Drug PK
      tabItem("pk",
        fluidRow(
          box(title="Cisplatin PK (central compartment)", width=6, status="warning",
            plotlyOutput("plot_cis_pk", height=280)
          ),
          box(title="Paclitaxel PK (central compartment, R/M)", width=6, status="warning",
            plotlyOutput("plot_pac_pk", height=280)
          )
        ),
        fluidRow(
          box(title="Pembrolizumab PK", width=6, status="danger",
            plotlyOutput("plot_pembro_pk", height=280)
          ),
          box(title="Tisotumab vedotin PK + free MMAE", width=6, status="danger",
            plotlyOutput("plot_tv_pk", height=280)
          )
        ),
        fluidRow(
          box(title="Bevacizumab PK + free VEGF-A", width=12, status="primary",
            plotlyOutput("plot_bev_vegf", height=300)
          )
        )
      ),

      ## Tab 3: PD Key Indicators
      tabItem("pd",
        fluidRow(
          box(title="Platinum-DNA adduct (relative)", width=6, status="warning",
            plotlyOutput("plot_ptdna", height=300)
          ),
          box(title="Cumulative radiation damage (LQ model)", width=6, status="danger",
            plotlyOutput("plot_rt", height=300)
          )
        ),
        fluidRow(
          box(title="CD8+ T cells (tumour-infiltrating CTL)", width=6, status="success",
            plotlyOutput("plot_cd8t", height=300)
          ),
          box(title="Tumour PD-L1 expression (CPS-like)", width=6, status="info",
            plotlyOutput("plot_pdl1", height=300)
          )
        ),
        fluidRow(
          box(title="PD indicator summary (final time point)", width=12, status="info",
            tableOutput("pd_summary")
          )
        )
      ),

      ## Tab 4: Clinical Endpoints
      tabItem("clinical",
        fluidRow(
          valueBoxOutput("box_best_resp"), valueBoxOutput("box_sccag_nadir"),
          valueBoxOutput("box_pfs")
        ),
        fluidRow(
          box(title="Tumour volume over time (cm³)", width=8, status="danger",
            plotlyOutput("plot_tv", height=380)
          ),
          box(title="RECIST response classification", width=4, status="warning",
            tableOutput("recist_table")
          )
        ),
        fluidRow(
          box(title="Tumour volume change rate (% vs baseline)", width=12, status="info",
            plotlyOutput("plot_tv_pct", height=250)
          )
        )
      ),

      ## Tab 5: Scenario Comparison
      tabItem("scenario",
        fluidRow(
          box(title="Six treatment scenarios — tumour volume comparison", width=12, status="danger",
            plotlyOutput("plot_scenario_tv", height=380)
          )
        ),
        fluidRow(
          box(title="Six treatment scenarios — SCC-Ag comparison", width=12, status="warning",
            plotlyOutput("plot_scenario_sccag", height=300)
          )
        ),
        fluidRow(
          box(title="Scenario summary table", width=12, status="info",
            DTOutput("scenario_table")
          )
        )
      ),

      ## Tab 6: Biomarkers
      tabItem("biomarker",
        fluidRow(
          box(title="Composite biomarker panel (time-concentration)", width=12, status="primary",
            plotlyOutput("plot_bm_panel", height=500)
          )
        ),
        fluidRow(
          box(title="Cervical cancer treatment algorithm (by stage)", width=6, status="danger",
            HTML('
<div style="font-family:monospace;line-height:1.8;font-size:13px">
<b>Cervical cancer treatment algorithm</b><br>
└─ <b>Stage classification (FIGO 2018)</b><br>
   ├─ <b>IB1-IB2 (localised)</b><br>
   │  └─ Radical hysterectomy ± adjuvant radiation<br>
   ├─ <b>IB3-IVA (locally advanced)</b><br>
   │  ├─ Cisplatin concurrent chemoradiotherapy (CCRT)<br>
   │  └─ High-risk group: +Pembrolizumab (KEYNOTE-A18)<br>
   └─ <b>IVB/recurrent/metastatic</b><br>
      ├─ 1st line: platinum+Paclitaxel+Bev±Pembro (GOG-240/KEYNOTE-826)<br>
      └─ 2nd line or later: Tisotumab vedotin (innovaTV 204/301)
</div>'
            )
          ),
          box(title="Key clinical trial reference values", width=6, status="warning",
            tableOutput("trial_ref_table")
          )
        )
      )
    )
  )
)

## ── Server ─────────────────────────────────────────────────────────────────
server <- function(input, output, session) {

  ## Reactive: build dosing events and run simulation
  sim_result <- eventReactive(input$run_sim, {
    TV0_v    <- input$TV0
    SCCAg0_v <- input$SCCAg_0
    cps_v    <- as.numeric(input$CPS_high)
    dur      <- input$sim_days
    rt_f     <- if(input$add_RT) 1 else 0
    ici_f    <- if(input$add_pembro) 1 else 0
    n_cis    <- input$n_cis

    ev_list <- list()

    if(input$add_RT) {
      cis_times <- seq(0, (n_cis-1)*7, by=7)
      ev_list <- c(ev_list, list(data.frame(time=cis_times, cmt="CIS_C1", amt=68/15, evid=1)))
    } else {
      ## Recurrent/metastatic interval chemo: cisplatin + paclitaxel q3w
      cis_times <- seq(0, 5*21, by=21)
      ev_list <- c(ev_list, list(data.frame(time=cis_times, cmt="CIS_C1", amt=68/15, evid=1)))
      ev_list <- c(ev_list, list(data.frame(time=cis_times, cmt="PAC_C1", amt=300000/6.5, evid=1)))
    }

    if(input$add_bev) {
      bev_times <- seq(0, dur, by=21)
      ev_list <- c(ev_list, list(data.frame(time=bev_times, cmt="BEV_C1", amt=1050/2.91, evid=1)))
    }
    if(input$add_pembro) {
      pembro_times <- seq(0, dur, by=21)
      ev_list <- c(ev_list, list(data.frame(time=pembro_times, cmt="PEMBRO_C1", amt=200/3.34, evid=1)))
    }
    if(input$add_tv) {
      tv_times <- seq(0, dur, by=21)
      ev_list <- c(ev_list, list(data.frame(time=tv_times, cmt="TV_ADC_C1", amt=140/3.35, evid=1)))
    }

    ev_df <- bind_rows(ev_list) %>% arrange(time)

    mod_run <- mod_base %>%
      param(TV0=TV0_v, SCCAg_0=SCCAg0_v, CPS_high=cps_v,
            ICI_flag=ici_f, RT_flag=rt_f)

    out <- mrgsim_df(mod_run, data=as_data_set(as.ev(ev_df)),
                     end=dur, delta=1)
    out
  }, ignoreNULL=FALSE)

  ## Six-scenario simulation (fixed for comparison tab)
  all_scenarios <- reactive({
    scenarios <- list(
      list(name="S1: Untreated", rt=0, ici=0, bev=FALSE, tv=FALSE, chemo_rm=FALSE),
      list(name="S2: Cisplatin CCRT", rt=1, ici=0, bev=FALSE, tv=FALSE, chemo_rm=FALSE),
      list(name="S3: CCRT+Pembrolizumab", rt=1, ici=1, bev=FALSE, tv=FALSE, chemo_rm=FALSE),
      list(name="S4: Chemo+Bev (GOG-240)", rt=0, ici=0, bev=TRUE, tv=FALSE, chemo_rm=TRUE),
      list(name="S5: Tisotumab vedotin", rt=0, ici=0, bev=FALSE, tv=TRUE, chemo_rm=FALSE),
      list(name="S6: Chemo+Bev+Pembro (KEYNOTE-826)", rt=0, ici=1, bev=TRUE, tv=FALSE, chemo_rm=TRUE)
    )

    lapply(scenarios, function(s) {
      ev_list <- list()
      if(s$rt == 1) {
        ev_list <- c(ev_list, list(data.frame(time=seq(0,5*7,by=7), cmt="CIS_C1", amt=68/15, evid=1)))
      } else if(s$chemo_rm) {
        ev_list <- c(ev_list, list(data.frame(time=seq(0,5*21,by=21), cmt="CIS_C1", amt=68/15, evid=1)))
        ev_list <- c(ev_list, list(data.frame(time=seq(0,5*21,by=21), cmt="PAC_C1", amt=300000/6.5, evid=1)))
      }
      if(s$bev) {
        ev_list <- c(ev_list, list(data.frame(time=seq(0,730,by=21), cmt="BEV_C1", amt=1050/2.91, evid=1)))
      }
      if(s$ici == 1) {
        ev_list <- c(ev_list, list(data.frame(time=seq(0,730,by=21), cmt="PEMBRO_C1", amt=200/3.34, evid=1)))
      }
      if(s$tv) {
        ev_list <- c(ev_list, list(data.frame(time=seq(0,730,by=21), cmt="TV_ADC_C1", amt=140/3.35, evid=1)))
      }
      ev_df <- bind_rows(ev_list) %>% arrange(time)
      mod_run <- mod_base %>% param(RT_flag=s$rt, ICI_flag=s$ici, CPS_high=1)
      out <- mrgsim_df(mod_run, data=as_data_set(as.ev(ev_df)), end=730, delta=1)
      out$Scenario <- s$name
      out
    }) %>% bind_rows()
  })

  ## ── ValueBoxes ───────────────────────────────────────────────
  output$box_stage <- renderValueBox({
    label <- switch(input$FIGO, early="FIGO IB1-IB2", lacc="FIGO IIB-IVA", rm="FIGO IVB/recurrent")
    valueBox(label, "Stage (FIGO 2018)", icon=icon("hospital"), color="red")
  })
  output$box_cps <- renderValueBox({
    label <- if(input$CPS_high=="1") "CPS ≥ 1 (positive)" else "CPS < 1 (negative)"
    valueBox(label, "PD-L1 CPS status", icon=icon("dna"), color="purple")
  })
  output$box_sccag0 <- renderValueBox({
    valueBox(paste0(input$SCCAg_0, " ng/mL"), "Initial SCC-Ag", icon=icon("vial"), color="orange")
  })
  output$box_best_resp <- renderValueBox({
    df <- sim_result()
    best <- min(df$TumorVol, na.rm=TRUE)
    tv0  <- input$TV0
    pct  <- round((best - tv0)/tv0*100, 1)
    valueBox(paste0(pct, "%"), "Best tumour response (RECIST)", icon=icon("chart-line"), color="red")
  })
  output$box_sccag_nadir <- renderValueBox({
    df <- sim_result()
    nadir <- round(min(df$SCCAg_lvl, na.rm=TRUE), 2)
    valueBox(paste0(nadir, " ng/mL"), "SCC-Ag nadir", icon=icon("flask"), color="purple")
  })
  output$box_pfs <- renderValueBox({
    df <- sim_result()
    tv_thresh <- input$TV0 * 2
    pfs_d <- df %>% filter(TumorVol > tv_thresh) %>% pull(time) %>% min()
    label <- if(is.infinite(pfs_d) || is.na(pfs_d))
      paste0(">", input$sim_days, " days")
    else
      paste0(round(pfs_d), " days (", round(pfs_d/30.4, 1), " mo)")
    valueBox(label, "PFS estimate (tumour-doubling time point)", icon=icon("calendar"), color="green")
  })

  ## ── Tab 1: Patient Summary ────────────────────────────────────
  output$patient_summary <- renderTable({
    data.frame(
      Item = c("Initial tumour size","Stage","PD-L1 CPS","Initial SCC-Ag","CCRT given","Number of cisplatin doses"),
      Value = c(paste(input$TV0,"cm³"),
               switch(input$FIGO, early="IB1-IB2", lacc="IIB-IVA", rm="IVB/recurrent"),
               ifelse(input$CPS_high=="1","≥1 (positive)","<1 (negative)"),
               paste(input$SCCAg_0,"ng/mL"),
               ifelse(input$add_RT,"Yes","No (R/M intermittent chemotherapy)"),
               paste(input$n_cis,"doses"))
    )
  })

  output$eligibility_table <- renderTable({
    cps <- input$CPS_high == "1"
    data.frame(
      Treatment = c("Pembrolizumab + CCRT (KEYNOTE-A18)", "Pembrolizumab 1st-line R/M (KEYNOTE-826)",
               "Bevacizumab (GOG-240)", "Tisotumab vedotin (innovaTV 301, 2L+)"),
      Eligibility = c(
        "✅ High-risk locally advanced (any stage)",
        ifelse(cps, "✅ Eligible (CPS≥1)", "⚠️ Limited benefit (CPS<1)"),
        "✅ All R/M patients (absent contraindications)",
        "✅ Eligible on progression after platinum therapy"
      )
    )
  })

  ## ── Tab 2: PK Plots ────────────────────────────────────────────
  pk_plot <- function(df, col, title, ylab, xlim=NULL) {
    p <- ggplot(df, aes(x=time, y=.data[[col]])) +
      geom_line(color="#C0392B") +
      labs(title=title, x="Day", y=ylab) +
      theme_bw(base_size=10)
    if(!is.null(xlim)) p <- p + coord_cartesian(xlim=xlim)
    ggplotly(p)
  }

  output$plot_cis_pk <- renderPlotly({
    df <- sim_result()
    pk_plot(df, "CIS_Conc", "Cisplatin PK", "Concentration (µg/mL)", xlim=c(0,120))
  })
  output$plot_pac_pk <- renderPlotly({
    df <- sim_result()
    pk_plot(df, "PAC_Conc", "Paclitaxel PK (R/M)", "Concentration (ng/mL)", xlim=c(0,250))
  })
  output$plot_pembro_pk <- renderPlotly({
    df <- sim_result()
    pk_plot(df, "PEMBRO_Conc", "Pembrolizumab PK", "Concentration (µg/mL)")
  })
  output$plot_tv_pk <- renderPlotly({
    df <- sim_result()
    p <- ggplot(df) +
      geom_line(aes(x=time, y=TVADC_Conc, color="TV-ADC (µg/mL)")) +
      geom_line(aes(x=time, y=MMAE_lvl, color="Free MMAE (relative)")) +
      labs(title="Tisotumab vedotin PK & MMAE release", x="Day", y="Concentration") +
      scale_color_manual(values=c("TV-ADC (µg/mL)"="#5E35B1","Free MMAE (relative)"="#FF7043")) +
      theme_bw(base_size=10) + theme(legend.position="bottom")
    ggplotly(p)
  })
  output$plot_bev_vegf <- renderPlotly({
    df <- sim_result()
    p <- ggplot(df) +
      geom_line(aes(x=time, y=BEV_Conc*10, color="Bevacizumab (×10 mg/L)")) +
      geom_line(aes(x=time, y=VEGF_free*1000, color="Free VEGF-A (×1000 ng/mL)")) +
      labs(title="Bevacizumab PK & VEGF suppression", x="Day", y="Concentration") +
      scale_color_manual(values=c("Bevacizumab (×10 mg/L)"="#1A237E","Free VEGF-A (×1000 ng/mL)"="#FF8A65")) +
      theme_bw(base_size=10) + theme(legend.position="bottom")
    ggplotly(p)
  })

  ## ── Tab 3: PD Key Indicators ───────────────────────────────────
  output$plot_ptdna <- renderPlotly({
    df <- sim_result()
    p <- ggplot(df, aes(x=time, y=PtDNA_rel)) +
      geom_line(color="#FF6F00") +
      labs(title="Platinum-DNA adduct (relative)", x="Day", y="Adduct (0–1)") +
      theme_bw(base_size=10)
    ggplotly(p)
  })
  output$plot_rt <- renderPlotly({
    df <- sim_result()
    p <- ggplot(df, aes(x=time, y=RT_damage)) +
      geom_line(color="#607D8B") +
      labs(title="Cumulative radiation damage (LQ model)", x="Day", y="Damage (relative)") +
      coord_cartesian(xlim=c(0,150)) +
      theme_bw(base_size=10)
    ggplotly(p)
  })
  output$plot_cd8t <- renderPlotly({
    df <- sim_result()
    p <- ggplot(df, aes(x=time, y=CD8T_rel)) +
      geom_line(color="#2E7D32") +
      geom_hline(yintercept=1.0, linetype="dashed", color="gray50") +
      labs(title="CD8+ T cells (intratumoural CTL)", x="Day", y="Relative level (baseline=1)") +
      theme_bw(base_size=10)
    ggplotly(p)
  })
  output$plot_pdl1 <- renderPlotly({
    df <- sim_result()
    p <- ggplot(df, aes(x=time, y=PDL1_rel)) +
      geom_line(color="#8E24AA") +
      labs(title="Tumour PD-L1 expression (CPS-like indicator)", x="Day", y="Relative expression (baseline=1)") +
      theme_bw(base_size=10)
    ggplotly(p)
  })
  output$pd_summary <- renderTable({
    df <- sim_result()
    last <- tail(df, 1)
    data.frame(
      Metric = c("Platinum-DNA adduct (relative)","Cumulative RT damage (relative)","CD8+ T cells (relative)","PD-L1 expression (relative)","Free VEGF-A (ng/mL)","Tumour volume (cm³)"),
      FinalValue = c(round(last$PtDNA_rel,3), round(last$RT_damage,3),
                  round(last$CD8T_rel,3), round(last$PDL1_rel,3),
                  round(last$VEGF_free,3), round(last$TumorVol,1))
    )
  })

  ## ── Tab 4: Clinical Endpoints ─────────────────────────────────
  output$plot_tv <- renderPlotly({
    df <- sim_result()
    thresh <- input$TV0 * 2
    p <- ggplot(df, aes(x=time, y=TumorVol)) +
      geom_line(color="#C0392B", size=1) +
      geom_hline(yintercept=thresh, linetype="dashed", color="gray40") +
      annotate("text", x=max(df$time)*0.7, y=thresh*1.05,
               label=paste("PD threshold (2×baseline,", round(thresh,0), "cm³)"), size=3, color="gray40") +
      labs(title="Tumour volume (cm³)", x="Day", y="Tumour volume (cm³)") +
      theme_bw(base_size=10)
    ggplotly(p)
  })
  output$plot_tv_pct <- renderPlotly({
    df <- sim_result()
    tv0 <- input$TV0
    df$pct_chg <- (df$TumorVol - tv0) / tv0 * 100
    p <- ggplot(df, aes(x=time, y=pct_chg)) +
      geom_line(color="#E65100") +
      geom_hline(yintercept=c(-30, 0, 20), linetype="dashed",
                 color=c("#27AE60","black","#C0392B")) +
      labs(title="Tumour volume change rate (% vs baseline)",
           x="Day", y="Change rate (%)") +
      theme_bw(base_size=10)
    ggplotly(p)
  })
  output$recist_table <- renderTable({
    df <- sim_result()
    tv0 <- input$TV0
    min_tv <- min(df$TumorVol)
    pct <- (min_tv - tv0) / tv0 * 100
    cat_r <- if(pct < -100) "CR (complete response)" else if(pct < -30) "PR (partial response)" else if(pct < 20) "SD (stable disease)" else "PD (progressive disease)"
    data.frame(
      Category   = c("Best response %", "RECIST category", "SCC-Ag nadir"),
      Value     = c(paste0(round(pct,1), "%"), cat_r, paste0(round(min(df$SCCAg_lvl),2), " ng/mL"))
    )
  })

  ## ── Tab 5: Scenario Comparison ────────────────────────────────
  output$plot_scenario_tv <- renderPlotly({
    df <- all_scenarios()
    p <- ggplot(df, aes(x=time, y=TumorVol, color=Scenario)) +
      geom_line(size=0.8) +
      geom_hline(yintercept=80, linetype="dashed", color="gray50") +
      labs(title="Six scenarios — tumour volume comparison (2 years)",
           x="Day", y="Tumour volume (cm³)") +
      scale_color_brewer(palette="Set1") +
      theme_bw(base_size=10) + theme(legend.position="bottom")
    ggplotly(p) %>% layout(legend=list(orientation="h", y=-0.3))
  })
  output$plot_scenario_sccag <- renderPlotly({
    df <- all_scenarios()
    p <- ggplot(df, aes(x=time, y=SCCAg_lvl, color=Scenario)) +
      geom_line(size=0.8) +
      geom_hline(yintercept=1.5, linetype="dashed", color="darkgreen") +
      scale_y_log10() +
      labs(title="Six scenarios — SCC-Ag (log scale)",
           x="Day", y="SCC-Ag (ng/mL)") +
      scale_color_brewer(palette="Set1") +
      theme_bw(base_size=10) + theme(legend.position="bottom")
    ggplotly(p) %>% layout(legend=list(orientation="h", y=-0.3))
  })
  output$scenario_table <- renderDT({
    df <- all_scenarios()
    df %>% group_by(Scenario) %>% summarise(
      MinTumour_cm3 = round(min(TumorVol, na.rm=TRUE), 1),
      MinSCCAg_ng_mL = round(min(SCCAg_lvl, na.rm=TRUE), 2),
      MaxRTDamage = round(max(RT_damage, na.rm=TRUE), 3),
      FinalCD8T = round(last(CD8T_rel), 3)
    ) %>% datatable(options=list(pageLength=6, dom="t"), rownames=FALSE)
  })

  ## ── Tab 6: Biomarker Panel ────────────────────────────────────
  output$plot_bm_panel <- renderPlotly({
    df <- sim_result()
    df_long <- df %>%
      select(time, SCCAg_lvl, TumorVol, PtDNA_rel, HPV_rel, CD8T_rel, VEGF_free) %>%
      pivot_longer(-time, names_to="Biomarker", values_to="Value") %>%
      mutate(Biomarker = recode(Biomarker,
        SCCAg_lvl="SCC-Ag (ng/mL)", TumorVol="Tumour volume (cm³)",
        PtDNA_rel="Pt-DNA adduct", HPV_rel="HPV viral load",
        CD8T_rel="CD8+ T cells", VEGF_free="Free VEGF-A (ng/mL)"))
    p <- ggplot(df_long, aes(x=time, y=Value)) +
      geom_line(aes(color=Biomarker)) +
      facet_wrap(~Biomarker, scales="free_y", ncol=3) +
      labs(title="Composite biomarker panel", x="Day", y="Value") +
      theme_bw(base_size=9) + theme(legend.position="none")
    ggplotly(p)
  })

  output$trial_ref_table <- renderTable({
    data.frame(
      Trial = c("RTOG-90-01 (1999/2004)","GOG-240 (2014 NEJM)","KEYNOTE-A18 (2024 Lancet)","KEYNOTE-826 (2021 NEJM)","innovaTV 301 (2024)"),
      Treatment = c("Cisplatin+5FU+RT vs extended-field RT","Chemo±Bevacizumab (R/M)","Pembrolizumab+CCRT vs CCRT (locally advanced)","Chemo±Bev±Pembro 1st-line (R/M)","Tisotumab vedotin vs chemotherapy (2L+)"),
      KeyResult = c("OS improvement (chemoradiotherapy superior)","mOS 16.8 vs 13.3mo (HR 0.71)","PFS HR 0.70 (high-risk group)","mOS improvement, greatest benefit in CPS≥1 group","mOS 11.5 vs 9.5mo (HR 0.70)"),
      PrimaryPopulation = c("Locally advanced","Recurrent/metastatic","Locally advanced, high-risk","1st-line recurrent/metastatic","2nd-line-plus recurrent/metastatic")
    )
  })
}

shinyApp(ui, server)
