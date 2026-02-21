import { useEffect, useRef, useState } from 'react'
import {
  Dna,
  FlaskConical,
  Brain,
  Filter,
  Activity,
  CheckCircle2,
  ArrowRight,
  ChevronDown,
  Pill,
  Leaf,
  Sparkles,
  Target,
  Heart,
  Zap,
  TrendingUp,
  Shield,
  Flame,
  HelpCircle,
  X,
  TrendingDown
} from 'lucide-react'

// Method explanation popup component
function MethodTooltip({ title, description, isOpen, onClose }) {
  if (!isOpen) return null

  return (
    <div className="fixed inset-0 z-50 flex items-center justify-center p-4 bg-black/50" onClick={onClose}>
      <div className="bg-white rounded-2xl p-6 max-w-md shadow-2xl" onClick={e => e.stopPropagation()}>
        <div className="flex justify-between items-start mb-4">
          <h3 className="font-bold text-xl text-gray-900">{title}</h3>
          <button onClick={onClose} className="text-gray-400 hover:text-gray-600">
            <X className="w-5 h-5" />
          </button>
        </div>
        <p className="text-gray-600">{description}</p>
      </div>
    </div>
  )
}

function App() {
  const observerRef = useRef(null)
  const [openTooltip, setOpenTooltip] = useState(null)

  useEffect(() => {
    observerRef.current = new IntersectionObserver(
      (entries) => {
        entries.forEach((entry) => {
          if (entry.isIntersecting) {
            entry.target.classList.add('animate-fade-in-up')
            entry.target.style.opacity = '1'
          }
        })
      },
      { threshold: 0.1 }
    )

    document.querySelectorAll('.fade-in').forEach((el) => {
      el.style.opacity = '0'
      observerRef.current.observe(el)
    })

    return () => observerRef.current?.disconnect()
  }, [])

  const scrollToSection = (id) => {
    document.getElementById(id)?.scrollIntoView({ behavior: 'smooth' })
  }

  const methodExplanations = {
    gnina: {
      title: "GNINA Docking",
      description: "A deep learning-enhanced molecular docking program that predicts how small molecules bind to proteins. Uses convolutional neural networks to score protein-ligand poses, achieving better accuracy than traditional scoring functions."
    },
    ml: {
      title: "Machine Learning (NOD2-Scout)",
      description: "A custom XGBoost classifier trained on NOD2-specific features. Filters docking results by learning patterns that distinguish true binders from false positives. Achieved 0.85-0.93 AUC-ROC on 5-fold cross-validation with scaffold split."
    },
    admet: {
      title: "ADMET Filtering",
      description: "Absorption, Distribution, Metabolism, Excretion, Toxicity screening. Applies drug-likeness rules (Lipinski's Rule of 5, PAINS filters) to ensure compounds can actually work as drugs in the human body."
    },
    md: {
      title: "Molecular Dynamics",
      description: "Physics-based simulation of protein-ligand complexes over time. Uses force fields to model atomic interactions. Confirms that compounds remain bound stably in the pocket rather than drifting away."
    },
    fep: {
      title: "Free Energy Perturbation",
      description: "Gold-standard method for calculating binding affinity differences. Uses thermodynamic cycles to compute ΔΔG between wild-type and mutant proteins. Computationally expensive but highly accurate."
    }
  }

  return (
    <div className="min-h-screen bg-white text-gray-900">
      {/* Method Tooltips */}
      {Object.entries(methodExplanations).map(([key, { title, description }]) => (
        <MethodTooltip
          key={key}
          title={title}
          description={description}
          isOpen={openTooltip === key}
          onClose={() => setOpenTooltip(null)}
        />
      ))}

      {/* Hero Section */}
      <section className="min-h-screen flex flex-col justify-center items-center px-6 bg-gradient-to-b from-slate-50 to-white relative">
        <div className="max-w-4xl mx-auto text-center">
          <div className="inline-flex items-center gap-2 px-4 py-2 bg-teal-50 text-teal-700 rounded-full text-sm font-medium mb-8">
            <Dna className="w-4 h-4" />
            Computational Drug Discovery
          </div>

          <h1 className="text-4xl md:text-6xl font-bold text-gray-900 mb-6 leading-tight">
            Precision Medicine, N=1
          </h1>

          <p className="text-xl md:text-2xl text-gray-600 mb-4">
            An FEP-Validated Deep Learning Pipeline to Rescue the NOD2 R702W Crohn's Variant
          </p>

          <p className="text-lg text-teal-600 font-medium mb-12 max-w-2xl mx-auto">
            A computational drug discovery pipeline for Crohn's disease patients with NOD2 mutations
          </p>

          <button
            onClick={() => scrollToSection('problem')}
            className="inline-flex items-center gap-2 px-8 py-4 bg-teal-600 text-white rounded-full font-semibold hover:bg-teal-700 transition-all hover:scale-105 shadow-lg shadow-teal-200"
          >
            Explore the Research
            <ChevronDown className="w-5 h-5" />
          </button>
        </div>

        <div className="absolute bottom-8 left-1/2 -translate-x-1/2 animate-bounce">
          <ChevronDown className="w-6 h-6 text-gray-400" />
        </div>
      </section>

      {/* What is Crohn's & NOD2 */}
      <section id="problem" className="py-24 px-6 bg-white">
        <div className="max-w-6xl mx-auto">
          {/* Explainer Cards */}
          <div className="fade-in grid md:grid-cols-2 gap-8 mb-20">
            {/* What is Crohn's */}
            <div className="bg-gradient-to-br from-orange-50 to-red-50 p-8 rounded-3xl border border-orange-100">
              <h3 className="text-2xl font-bold mb-4 text-gray-900">What is Crohn's Disease?</h3>
              <p className="text-gray-700 mb-4">
                Crohn's disease is a chronic inflammatory bowel disease (IBD) that causes inflammation in the digestive tract. It affects <strong>3 million Americans</strong> and has no cure.
              </p>
              <ul className="space-y-2 text-gray-600">
                <li className="flex items-start gap-2">
                  <span className="text-orange-500 mt-1">•</span>
                  <span>Causes severe abdominal pain, diarrhea, fatigue, and malnutrition</span>
                </li>
                <li className="flex items-start gap-2">
                  <span className="text-orange-500 mt-1">•</span>
                  <span>Often requires surgery - 70% of patients need it within their lifetime</span>
                </li>
                <li className="flex items-start gap-2">
                  <span className="text-orange-500 mt-1">•</span>
                  <span>Current treatments suppress immunity but don't target the root cause</span>
                </li>
              </ul>
            </div>

            {/* What is NOD2 */}
            <div className="bg-gradient-to-br from-blue-50 to-indigo-50 p-8 rounded-3xl border border-blue-100">
              <h3 className="text-2xl font-bold mb-4 text-gray-900">What is NOD2?</h3>
              <p className="text-gray-700 mb-4">
                NOD2 is an immune sensor protein that detects bacterial invaders in the gut. When it's mutated, the immune system malfunctions.
              </p>
              <ul className="space-y-2 text-gray-600">
                <li className="flex items-start gap-2">
                  <span className="text-blue-500 mt-1">•</span>
                  <span>Mutations in NOD2 are the <strong>#1 genetic risk factor</strong> for Crohn's</span>
                </li>
                <li className="flex items-start gap-2">
                  <span className="text-blue-500 mt-1">•</span>
                  <span>R702W is the most common mutation (found in ~15% of Crohn's patients)</span>
                </li>
                <li className="flex items-start gap-2">
                  <span className="text-blue-500 mt-1">•</span>
                  <span>No drugs currently target NOD2 directly - that's the gap I'm filling</span>
                </li>
              </ul>
            </div>
          </div>

          <div className="fade-in text-center mb-16">
            <h2 className="text-3xl md:text-4xl font-bold mb-4">The Problem</h2>
            <p className="text-gray-600 text-lg max-w-2xl mx-auto">
              Crohn's disease affects millions, but genetic mutations make treatment complex
            </p>
          </div>

          <div className="grid md:grid-cols-3 gap-8 mb-16">
            <div className="fade-in bg-gradient-to-br from-red-50 to-orange-50 p-8 rounded-2xl border border-red-100">
              <div className="w-14 h-14 bg-red-100 rounded-xl flex items-center justify-center mb-6">
                <Dna className="w-7 h-7 text-red-600" />
              </div>
              <h3 className="text-xl font-bold mb-3">NOD2 Mutations</h3>
              <p className="text-gray-600 mb-4">
                The #1 genetic risk factor for Crohn's disease
              </p>
              <div className="text-5xl font-bold text-red-600">40×</div>
              <p className="text-sm text-gray-500 mt-1">increased risk</p>
            </div>

            <div className="fade-in bg-gradient-to-br from-amber-50 to-yellow-50 p-8 rounded-2xl border border-amber-100">
              <div className="w-14 h-14 bg-amber-100 rounded-xl flex items-center justify-center mb-6">
                <Target className="w-7 h-7 text-amber-600" />
              </div>
              <h3 className="text-xl font-bold mb-3">No Targeted Drugs</h3>
              <p className="text-gray-600 mb-4">
                No drugs specifically target the NOD2 protein
              </p>
              <div className="text-5xl font-bold text-amber-600">0</div>
              <p className="text-sm text-gray-500 mt-1">approved NOD2 drugs</p>
            </div>

            <div className="fade-in bg-gradient-to-br from-teal-50 to-cyan-50 p-8 rounded-2xl border border-teal-100">
              <div className="w-14 h-14 bg-teal-100 rounded-xl flex items-center justify-center mb-6">
                <Heart className="w-7 h-7 text-teal-600" />
              </div>
              <h3 className="text-xl font-bold mb-3">Personal Mission</h3>
              <p className="text-gray-600 mb-4">
                I have Crohn's disease with the R702W mutation
              </p>
              <div className="text-5xl font-bold text-teal-600">N=1</div>
              <p className="text-sm text-gray-500 mt-1">precision medicine</p>
            </div>
          </div>

          {/* Gene to Disease Flow */}
          <div className="fade-in flex flex-col md:flex-row items-center justify-center gap-4 md:gap-8 bg-slate-50 p-8 rounded-2xl">
            <div className="text-center">
              <div className="w-20 h-20 bg-blue-100 rounded-full flex items-center justify-center mx-auto mb-3 shadow-lg shadow-blue-100">
                <Dna className="w-10 h-10 text-blue-600" />
              </div>
              <p className="font-bold text-lg">NOD2 Gene</p>
              <p className="text-sm text-gray-500">Mutated (R702W)</p>
            </div>
            <ArrowRight className="w-8 h-8 text-gray-300 rotate-90 md:rotate-0" />
            <div className="text-center">
              <div className="w-20 h-20 bg-purple-100 rounded-full flex items-center justify-center mx-auto mb-3 shadow-lg shadow-purple-100">
                <FlaskConical className="w-10 h-10 text-purple-600" />
              </div>
              <p className="font-bold text-lg">NOD2 Protein</p>
              <p className="text-sm text-gray-500">Dysfunctional</p>
            </div>
            <ArrowRight className="w-8 h-8 text-gray-300 rotate-90 md:rotate-0" />
            <div className="text-center">
              <div className="w-20 h-20 bg-red-100 rounded-full flex items-center justify-center mx-auto mb-3 shadow-lg shadow-red-100">
                <Flame className="w-10 h-10 text-red-600" />
              </div>
              <p className="font-bold text-lg">Crohn's Disease</p>
              <p className="text-sm text-gray-500">Chronic inflammation</p>
            </div>
          </div>
        </div>
      </section>

      {/* The Pipeline */}
      <section className="py-24 px-6 bg-gradient-to-b from-slate-50 to-white">
        <div className="max-w-6xl mx-auto">
          <div className="fade-in text-center mb-16">
            <h2 className="text-3xl md:text-4xl font-bold mb-4">Screening Pipeline</h2>
            <p className="text-gray-600 text-lg max-w-2xl mx-auto">
              Multi-stage computational funnel from library to validated hits
            </p>
          </div>

          {/* Improved Pipeline Steps with ? explanations */}
          <div className="fade-in mb-16">
            <div className="grid grid-cols-2 md:grid-cols-3 lg:grid-cols-6 gap-4">
              {[
                {
                  icon: FlaskConical,
                  label: '9,566',
                  sublabel: 'Compounds',
                  color: 'blue',
                  desc: '7,414 natural + 2,152 FDA',
                  tooltip: null
                },
                {
                  icon: Target,
                  label: 'GNINA',
                  sublabel: 'Docking',
                  color: 'indigo',
                  desc: '→ 500 hits',
                  tooltip: 'gnina'
                },
                {
                  icon: Brain,
                  label: 'ML',
                  sublabel: 'NOD2-Scout',
                  color: 'violet',
                  desc: '0.85-0.93 AUC',
                  tooltip: 'ml'
                },
                {
                  icon: Filter,
                  label: 'ADMET',
                  sublabel: 'Filter',
                  color: 'purple',
                  desc: '→ 8 candidates',
                  tooltip: 'admet'
                },
                {
                  icon: Activity,
                  label: '520 ns',
                  sublabel: 'MD Sims',
                  color: 'pink',
                  desc: 'Binding stability',
                  tooltip: 'md'
                },
                {
                  icon: Zap,
                  label: 'FEP',
                  sublabel: '120 windows',
                  color: 'teal',
                  desc: '→ 2 validated',
                  tooltip: 'fep'
                },
              ].map((step, i) => {
                const colorMap = {
                  blue: { bg: 'bg-blue-50', border: 'border-blue-200', icon: 'text-blue-600', badge: 'bg-blue-100' },
                  indigo: { bg: 'bg-indigo-50', border: 'border-indigo-200', icon: 'text-indigo-600', badge: 'bg-indigo-100' },
                  violet: { bg: 'bg-violet-50', border: 'border-violet-200', icon: 'text-violet-600', badge: 'bg-violet-100' },
                  purple: { bg: 'bg-purple-50', border: 'border-purple-200', icon: 'text-purple-600', badge: 'bg-purple-100' },
                  pink: { bg: 'bg-pink-50', border: 'border-pink-200', icon: 'text-pink-600', badge: 'bg-pink-100' },
                  teal: { bg: 'bg-teal-50', border: 'border-teal-200', icon: 'text-teal-600', badge: 'bg-teal-100' },
                }
                const colors = colorMap[step.color]

                return (
                  <div key={i} className="relative">
                    <div className={`${colors.bg} ${colors.border} border-2 rounded-2xl p-5 text-center h-full transition-all hover:shadow-lg hover:scale-105`}>
                      {step.tooltip && (
                        <button
                          onClick={() => setOpenTooltip(step.tooltip)}
                          className="absolute top-2 right-2 text-gray-400 hover:text-gray-600"
                        >
                          <HelpCircle className="w-5 h-5" />
                        </button>
                      )}
                      <div className={`w-12 h-12 ${colors.badge} rounded-xl flex items-center justify-center mx-auto mb-3`}>
                        <step.icon className={`w-6 h-6 ${colors.icon}`} />
                      </div>
                      <div className="font-bold text-2xl text-gray-900">{step.label}</div>
                      <div className="text-sm font-medium text-gray-600 mb-2">{step.sublabel}</div>
                      <div className="text-xs text-gray-500">{step.desc}</div>
                    </div>
                    {i < 5 && (
                      <div className="hidden lg:block absolute top-1/2 -right-3 transform -translate-y-1/2 z-10">
                        <ArrowRight className="w-5 h-5 text-gray-300" />
                      </div>
                    )}
                  </div>
                )
              })}
            </div>
          </div>

          {/* Screening Funnel Chart */}
          <div className="fade-in">
            <div className="bg-white p-6 rounded-3xl border border-slate-200 shadow-sm">
              <img
                src="/charts/screening_funnel.png"
                alt="Figure 3: Computational Screening Funnel showing compound reduction from 9,566 to 2 validated hits"
                className="w-full max-w-4xl mx-auto rounded-lg"
              />
            </div>
          </div>
        </div>
      </section>

      {/* Key Results with Python Charts */}
      <section className="py-24 px-6 bg-white">
        <div className="max-w-6xl mx-auto">
          <div className="fade-in text-center mb-16">
            <h2 className="text-3xl md:text-4xl font-bold mb-4">Key Results</h2>
            <p className="text-gray-600 text-lg max-w-2xl mx-auto">
              FEP-validated binding free energies for two drug candidates
            </p>
          </div>

          {/* Python-generated Charts */}
          <div className="fade-in grid md:grid-cols-2 gap-8 mb-12">
            <div className="bg-white p-6 rounded-3xl border border-slate-200 shadow-sm">
              <img
                src="/charts/binding_energy.png"
                alt="Figure 1: FEP-Calculated Binding Affinities showing Febuxostat and Bufadienolide binding energies with error bars"
                className="w-full rounded-lg"
              />
            </div>
            <div className="bg-white p-6 rounded-3xl border border-slate-200 shadow-sm">
              <img
                src="/charts/ddg_comparison.png"
                alt="Figure 2: Mutation Effect showing ΔΔG comparison between Febuxostat (+2.34) and Bufadienolide (-0.44)"
                className="w-full rounded-lg"
              />
            </div>
          </div>

          {/* WT vs MUT Chart */}
          <div className="fade-in mb-12">
            <div className="bg-white p-6 rounded-3xl border border-slate-200 shadow-sm">
              <img
                src="/charts/wt_vs_mut.png"
                alt="Figure 4: Wild-Type vs Mutant Binding comparison"
                className="w-full max-w-3xl mx-auto rounded-lg"
              />
            </div>
          </div>

          {/* Summary Cards */}
          <div className="grid md:grid-cols-2 gap-8">
            {/* Febuxostat Card */}
            <div className="fade-in bg-gradient-to-br from-blue-50 to-indigo-50 p-8 rounded-3xl border border-blue-100 relative overflow-hidden">
              <div className="absolute top-0 right-0 w-32 h-32 bg-blue-100 rounded-full -translate-y-1/2 translate-x-1/2 opacity-50" />

              <div className="relative">
                <div className="flex items-center gap-3 mb-6">
                  <div className="w-14 h-14 bg-blue-100 rounded-xl flex items-center justify-center">
                    <Pill className="w-7 h-7 text-blue-600" />
                  </div>
                  <div>
                    <h3 className="text-2xl font-bold">Febuxostat</h3>
                    <span className="text-sm text-blue-600 font-medium">FDA-Approved (Gout)</span>
                  </div>
                </div>

                <div className="space-y-3 mb-6">
                  <div className="flex justify-between items-center py-2 border-b border-blue-100">
                    <span className="text-gray-600">ΔG<sub>bind</sub> (WT)</span>
                    <span className="font-bold text-xl text-gray-900">-10.36 ± 0.18</span>
                  </div>
                  <div className="flex justify-between items-center py-2 border-b border-blue-100">
                    <span className="text-gray-600">ΔG<sub>bind</sub> (R702W)</span>
                    <span className="font-bold text-xl text-gray-900">-8.02 ± 0.19</span>
                  </div>
                  <div className="flex justify-between items-center py-2 border-b border-blue-100">
                    <span className="text-gray-600">ΔΔG</span>
                    <span className="font-bold text-2xl text-red-600">+2.34 kcal/mol</span>
                  </div>
                </div>

                <div className="p-4 bg-red-100/50 rounded-xl">
                  <p className="text-sm text-red-800 flex items-center gap-2">
                    <TrendingDown className="w-4 h-4" />
                    <span><strong>~50× weaker</strong> binding in mutant patients</span>
                  </p>
                </div>
              </div>
            </div>

            {/* Bufadienolide Card */}
            <div className="fade-in bg-gradient-to-br from-teal-50 to-emerald-50 p-8 rounded-3xl border border-teal-100 relative overflow-hidden">
              <div className="absolute top-0 right-0 w-32 h-32 bg-teal-100 rounded-full -translate-y-1/2 translate-x-1/2 opacity-50" />

              <div className="absolute top-4 right-4 bg-teal-600 text-white text-xs font-bold px-3 py-1 rounded-full flex items-center gap-1">
                <Sparkles className="w-3 h-3" />
                NOVEL DISCOVERY
              </div>

              <div className="relative">
                <div className="flex items-center gap-3 mb-6">
                  <div className="w-14 h-14 bg-teal-100 rounded-xl flex items-center justify-center">
                    <Leaf className="w-7 h-7 text-teal-600" />
                  </div>
                  <div>
                    <h3 className="text-2xl font-bold">Bufadienolide</h3>
                    <span className="text-sm text-teal-600 font-medium">Natural Product (CID 10120)</span>
                  </div>
                </div>

                <div className="space-y-3 mb-6">
                  <div className="flex justify-between items-center py-2 border-b border-teal-100">
                    <span className="text-gray-600">ΔG<sub>bind</sub> (WT)</span>
                    <span className="font-bold text-xl text-gray-900">-15.22 ± 0.26</span>
                  </div>
                  <div className="flex justify-between items-center py-2 border-b border-teal-100">
                    <span className="text-gray-600">ΔG<sub>bind</sub> (R702W)</span>
                    <span className="font-bold text-xl text-gray-900">-15.66 ± 0.26</span>
                  </div>
                  <div className="flex justify-between items-center py-2 border-b border-teal-100">
                    <span className="text-gray-600">ΔΔG</span>
                    <span className="font-bold text-2xl text-teal-600">-0.44 kcal/mol</span>
                  </div>
                </div>

                <div className="p-4 bg-teal-100/50 rounded-xl">
                  <p className="text-sm text-teal-800 flex items-center gap-2">
                    <Shield className="w-4 h-4" />
                    <span><strong>Mutation-resistant</strong> - works for ALL patients!</span>
                  </p>
                </div>
              </div>
            </div>
          </div>
        </div>
      </section>

      {/* Novel Discovery */}
      <section className="py-24 px-6 bg-gradient-to-b from-teal-600 to-teal-700 text-white">
        <div className="max-w-4xl mx-auto">
          <div className="fade-in text-center mb-12">
            <div className="inline-flex items-center gap-2 px-4 py-2 bg-white/20 rounded-full text-sm font-medium mb-6">
              <Sparkles className="w-4 h-4" />
              Novel Discovery
            </div>
            <h2 className="text-3xl md:text-4xl font-bold mb-6">First Report of Bufadienolide-NOD2 Binding</h2>
            <p className="text-xl text-teal-100">
              Nobody has ever connected this natural compound class to the NOD2 protein before
            </p>
          </div>

          <div className="fade-in grid md:grid-cols-3 gap-6">
            <div className="bg-white/10 backdrop-blur-sm p-6 rounded-2xl border border-white/20 text-center">
              <Leaf className="w-12 h-12 mx-auto mb-4 text-teal-200" />
              <h3 className="font-bold text-lg mb-2">Bufadienolides</h3>
              <p className="text-teal-100 text-sm">Natural compounds from toads and plants</p>
            </div>
            <div className="bg-white/10 backdrop-blur-sm p-6 rounded-2xl border border-white/20 text-center">
              <Flame className="w-12 h-12 mx-auto mb-4 text-orange-300" />
              <h3 className="font-bold text-lg mb-2">NF-κB Inhibitors</h3>
              <p className="text-teal-100 text-sm">Known to block inflammatory signaling</p>
            </div>
            <div className="bg-white/10 backdrop-blur-sm p-6 rounded-2xl border border-white/20 text-center">
              <Dna className="w-12 h-12 mx-auto mb-4 text-blue-300" />
              <h3 className="font-bold text-lg mb-2">NOD2 → NF-κB</h3>
              <p className="text-teal-100 text-sm">NOD2 signals through NF-κB pathway</p>
            </div>
          </div>

          <div className="fade-in mt-12 bg-white/10 backdrop-blur-sm p-8 rounded-2xl border border-white/20 text-center">
            <h3 className="text-2xl font-bold mb-4">Dual Mechanism of Action</h3>
            <p className="text-lg text-teal-100 max-w-2xl mx-auto">
              Bufadienolide may work through <strong>two pathways</strong>: directly binding NOD2 AND inhibiting downstream NF-κB inflammation. This could make it especially effective for Crohn's disease.
            </p>
          </div>
        </div>
      </section>

      {/* Methods Summary */}
      <section className="py-24 px-6 bg-slate-50">
        <div className="max-w-6xl mx-auto">
          <div className="fade-in text-center mb-16">
            <h2 className="text-3xl md:text-4xl font-bold mb-4">Methods</h2>
            <p className="text-gray-600">Click the <HelpCircle className="w-4 h-4 inline" /> icons above for detailed explanations</p>
          </div>

          <div className="fade-in grid md:grid-cols-2 lg:grid-cols-4 gap-6">
            {[
              { title: 'Virtual Screening', value: 'GNINA', desc: 'Deep learning docking', icon: Target, tooltip: 'gnina' },
              { title: 'Machine Learning', value: 'XGBoost', desc: 'Custom NOD2-Scout', icon: Brain, tooltip: 'ml' },
              { title: 'MD Simulation', value: '520 ns', desc: 'All-atom dynamics', icon: Activity, tooltip: 'md' },
              { title: 'Free Energy', value: 'FEP + MBAR', desc: 'Gold-standard ΔΔG', icon: Zap, tooltip: 'fep' },
            ].map((method, i) => (
              <div
                key={i}
                className="bg-white p-6 rounded-2xl border border-slate-200 shadow-sm cursor-pointer hover:shadow-md transition-all"
                onClick={() => setOpenTooltip(method.tooltip)}
              >
                <div className="flex items-start justify-between mb-4">
                  <method.icon className="w-8 h-8 text-teal-600" />
                  <HelpCircle className="w-5 h-5 text-gray-400" />
                </div>
                <h3 className="font-bold text-lg mb-1">{method.title}</h3>
                <div className="text-2xl font-bold text-gray-900 mb-2">{method.value}</div>
                <p className="text-sm text-gray-500">{method.desc}</p>
              </div>
            ))}
          </div>

          {/* Computational Stats */}
          <div className="fade-in mt-12">
            <div className="bg-white p-8 rounded-3xl border border-slate-200 shadow-sm">
              <h3 className="text-lg font-bold text-gray-800 text-center mb-8">Computational Scale</h3>
              <div className="grid md:grid-cols-3 gap-8">
                <div className="text-center">
                  <div className="text-5xl font-bold text-indigo-600 mb-2">520</div>
                  <div className="text-gray-600 font-medium">nanoseconds</div>
                  <div className="text-sm text-gray-500 mt-1">total MD simulation</div>
                </div>
                <div className="text-center">
                  <div className="text-5xl font-bold text-amber-600 mb-2">120</div>
                  <div className="text-gray-600 font-medium">FEP windows</div>
                  <div className="text-sm text-gray-500 mt-1">alchemical perturbations</div>
                </div>
                <div className="text-center">
                  <div className="text-5xl font-bold text-pink-600 mb-2">9,566</div>
                  <div className="text-gray-600 font-medium">compounds</div>
                  <div className="text-sm text-gray-500 mt-1">screened computationally</div>
                </div>
              </div>
            </div>
          </div>
        </div>
      </section>

      {/* Footer */}
      <footer className="py-12 px-6 bg-slate-900 text-white">
        <div className="max-w-6xl mx-auto text-center">
          <div className="flex items-center justify-center gap-2 mb-4">
            <Dna className="w-6 h-6 text-teal-400" />
            <span className="font-bold text-xl">NOD2 Drug Discovery</span>
          </div>
          <p className="text-slate-400 mb-2">Precision Medicine, N=1</p>
          <p className="text-slate-500 text-sm">
            Computational Drug Discovery for Crohn's Disease
          </p>
        </div>
      </footer>
    </div>
  )
}

export default App
