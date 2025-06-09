import React, { useState, useEffect, useRef } from 'react';
import { motion } from 'framer-motion';
import { useNavigate } from 'react-router-dom';
import {
  FaCompass,
  FaBalanceScale,
  FaFlask,
  FaUserSecret,
  FaTasks,
  FaChevronDown,
  FaChevronUp,
  FaDna,
  FaCubes as FaMolecule,
  FaFlask as FaLab,
  FaSearch,
  FaRobot,
  FaChartLine,
  FaDatabase,
  FaUsers,
  FaGraduationCap,
  FaLightbulb,
  FaAtom,
  FaCog,
  FaPlay,
  FaBrain,
  FaNetworkWired,
  FaMicrochip,
  FaFingerprint,
  FaCubes,
  FaInfinity,
  FaChartPie,
  FaBullseye,
  FaLock,
  FaRocket,
  FaGlobeAmericas,
  FaStopwatch,
  FaShieldAlt,
  FaCheckCircle,
  FaSyncAlt,
  FaEye,
  FaCodeBranch,
  FaMagic,
  FaUpload,
  FaArrowUp
} from 'react-icons/fa';
import GlassyContainer from '../components/glassy_container/gc';
import { fadeInUpVariantStatic, fadeInDownVariants, fadeInLeftVariants } from '../components/animations/framerAnim';
import './Home.css';

// New Infographic Components
const StatsCircle = ({ value, label, color, maxValue = 100 }) => {
  const [animatedValue, setAnimatedValue] = useState(0);
  const percentage = typeof value === 'string' ? 100 : (value / maxValue) * 100;

  useEffect(() => {
    const timer = setTimeout(() => {
      setAnimatedValue(percentage);
    }, 500);
    return () => clearTimeout(timer);
  }, [percentage]);

  return (
    <div className="stats-circle">
      <svg viewBox="0 0 100 100" className="circular-chart">
        <circle
          className="circle-bg"
          cx="50"
          cy="50"
          r="40"
        />
        <circle
          className="circle"
          cx="50"
          cy="50"
          r="40"
          style={{
            '--circle-color': color,
            strokeDasharray: `${animatedValue * 2.51}, 251.2`,
            transition: 'stroke-dasharray 1.5s ease-in-out'
          }}
        />
        <text
          x="50"
          y="45"
          className="circle-text-value"
          textAnchor="middle"
          dominantBaseline="middle"
        >
          {value}
        </text>
        <text
          x="50"
          y="60"
          className="circle-text-label"
          textAnchor="middle"
          dominantBaseline="middle"
        >
          {label}
        </text>
      </svg>
    </div>
  );
};

const ProgressBar = ({ label, percentage, color, icon }) => {
  const [animatedWidth, setAnimatedWidth] = useState(0);

  useEffect(() => {
    const timer = setTimeout(() => {
      setAnimatedWidth(percentage);
    }, 500);
    return () => clearTimeout(timer);
  }, [percentage]);

  return (
    <div className="progress-bar-container">
      <div className="progress-bar-header">
        <div className="progress-bar-icon" style={{ color: color }}>
          {icon}
        </div>
        <span className="progress-bar-label">{label}</span>
        <span className="progress-bar-percentage">{percentage}%</span>
      </div>
      <div className="progress-bar-track">
        <div
          className="progress-bar-fill"
          style={{
            width: `${animatedWidth}%`,
            backgroundColor: color,
            transition: 'width 1.5s ease-in-out'
          }}
        />
      </div>
    </div>
  );
};

const ProcessFlow = ({ steps }) => {
  return (
    <div className="process-flow">
      {steps.map((step, index) => (
        <motion.div
          key={index}
          className="process-step"
          initial={{ opacity: 0, y: 20 }}
          animate={{ opacity: 1, y: 0 }}
          transition={{ delay: index * 0.2 }}
        >
          <div className="process-step-number">{index + 1}</div>
          <div className="process-step-icon" style={{ color: step.color }}>
            {step.icon}
          </div>
          <div className="process-step-content">
            <h6 className="process-step-title">{step.title}</h6>
            <p className="process-step-description">{step.description}</p>
          </div>
          {index < steps.length - 1 && (
            <div className="process-connector">
              <FaChevronDown />
            </div>
          )}
        </motion.div>
      ))}
    </div>
  );
};

const NetworkGraph = ({ nodes }) => {
  return (
    <div className="network-graph">
      <svg viewBox="0 0 400 300" className="network-svg">
        {/* Connection lines */}
        <g className="connections">
          <line x1="200" y1="150" x2="100" y2="80" className="connection-line" />
          <line x1="200" y1="150" x2="300" y2="80" className="connection-line" />
          <line x1="200" y1="150" x2="100" y2="220" className="connection-line" />
          <line x1="200" y1="150" x2="300" y2="220" className="connection-line" />
          <line x1="100" y1="80" x2="300" y2="80" className="connection-line" />
          <line x1="100" y1="220" x2="300" y2="220" className="connection-line" />
        </g>

        {/* Central node */}
        <circle cx="200" cy="150" r="30" className="network-node central" />
        <text x="200" y="155" textAnchor="middle" className="node-text">AI Core</text>

        {/* Peripheral nodes */}
        <circle cx="100" cy="80" r="20" className="network-node" />
        <text x="100" y="85" textAnchor="middle" className="node-text-small">Drug</text>

        <circle cx="300" cy="80" r="20" className="network-node" />
        <text x="300" y="85" textAnchor="middle" className="node-text-small">Protein</text>

        <circle cx="100" cy="220" r="20" className="network-node" />
        <text x="100" y="225" textAnchor="middle" className="node-text-small">ADMET</text>

        <circle cx="300" cy="220" r="20" className="network-node" />
        <text x="300" y="225" textAnchor="middle" className="node-text-small">Optimize</text>
      </svg>
    </div>
  );
};

const Home = () => {
  const navigate = useNavigate();
  const [expandedSections, setExpandedSections] = useState(new Set(['overview']));
  const [activeFeature, setActiveFeature] = useState(null);

  // Timeline state and refs
  const [scrollProgress, setScrollProgress] = useState(0);
  const [activeSection, setActiveSection] = useState('hero');
  const timelineRef = useRef(null);
  const heroRef = useRef(null);
  const quickActionsRef = useRef(null);
  const platformRef = useRef(null);
  const ctaRef = useRef(null);

  // Scroll handler for timeline progress and active section
  useEffect(() => {
    const handleScroll = () => {
      const scrollTop = window.pageYOffset;
      const docHeight = document.documentElement.scrollHeight - window.innerHeight;
      const progress = (scrollTop / docHeight) * 100;
      setScrollProgress(progress);

      // Determine active section based on scroll position
      const sections = [
        { id: 'hero', ref: heroRef },
        { id: 'quick-actions', ref: quickActionsRef },
        { id: 'platform', ref: platformRef },
        { id: 'cta', ref: ctaRef }
      ];

      let currentActive = 'hero';
      sections.forEach(section => {
        if (section.ref.current) {
          const rect = section.ref.current.getBoundingClientRect();
          if (rect.top <= window.innerHeight / 2 && rect.bottom >= window.innerHeight / 2) {
            currentActive = section.id;
          }
        }
      });

      setActiveSection(currentActive);
    };

    window.addEventListener('scroll', handleScroll);
    handleScroll(); // Initial call

    return () => window.removeEventListener('scroll', handleScroll);
  }, []);

  const toggleSection = (sectionId) => {
    const newExpanded = new Set(expandedSections);
    if (newExpanded.has(sectionId)) {
      newExpanded.delete(sectionId);
    } else {
      newExpanded.add(sectionId);
    }
    setExpandedSections(newExpanded);
  };

  // Enhanced platform sections with infographics
  const platformSections = [
    {
      id: 'overview',
      title: 'Platform Overview',
      icon: <FaLightbulb />,
      color: '#22c55e',
      description: 'Comprehensive drug discovery and molecular analysis platform',
      content: {
        subtitle: 'Advanced AI-Powered Drug Discovery Toolkit',
        features: [
          'Integrated molecular analysis and visualization',
          'AI-driven optimization and generation algorithms',
          'Comprehensive database integration (ChEMBL, PubChem, UniProt)',
          'Real-time ADMET property prediction',
          'Protein structure analysis and prediction',
          'Interactive 3D molecular visualization'
        ],
        stats: [
          { label: 'Molecular Databases', value: '3+', icon: <FaDatabase /> },
          { label: 'Analysis Tools', value: '15+', icon: <FaCog /> },
          { label: 'AI Agents', value: '4', icon: <FaRobot /> },
          { label: 'Research Areas', value: '∞', icon: <FaGraduationCap /> }
        ],
        infographics: {
          circularStats: [
            { value: '95%', label: 'Accuracy', color: '#22c55e', maxValue: 100 },
            { value: '99.9%', label: 'Uptime', color: '#16a34a', maxValue: 100 },
            { value: '24/7', label: 'Support', color: '#15803d', maxValue: 100 }
          ],
          processSteps: [
            {
              title: 'Data Input',
              description: 'Submit molecular structures or targets',
              icon: <FaUpload />,
              color: '#22c55e'
            },
            {
              title: 'AI Processing',
              description: 'Advanced algorithms analyze your data',
              icon: <FaBrain />,
              color: '#16a34a'
            },
            {
              title: 'Results',
              description: 'Get comprehensive insights and predictions',
              icon: <FaChartLine />,
              color: '#15803d'
            }
          ]
        }
      }
    },
    {
      id: 'master-agent',
      title: 'Master Agent System',
      icon: <FaUserSecret />,
      color: '#22c55e',
      description: 'Intelligent AI orchestrator for complex drug discovery workflows',
      content: {
        subtitle: 'Your AI Research Assistant',
        features: [
          'Natural language query processing',
          'Multi-agent coordination and delegation',
          'Context-aware task routing',
          'Specialized sub-agents for different domains',
          'Real-time workflow management',
          'Interactive conversational interface'
        ],
        agents: [
          { name: 'DrugAgent', specialty: 'Database searches & exploration', icon: <FaSearch /> },
          { name: 'DrugOptimizationAgent', specialty: 'ADMET & molecular generation', icon: <FaFlask /> },
          { name: 'ProteinAgent', specialty: 'Protein analysis & structure', icon: <FaDna /> }
        ],
        infographics: {
          networkGraph: true,
          agentMetrics: [
            { label: 'Response Time', percentage: 98, color: '#22c55e', icon: <FaStopwatch /> },
            { label: 'Query Accuracy', percentage: 94, color: '#16a34a', icon: <FaBullseye /> },
            { label: 'Task Success', percentage: 96, color: '#15803d', icon: <FaCheckCircle /> }
          ]
        }
      }
    },
    {
      id: 'structure-analysis',
      title: 'Structure Analysis',
      icon: <FaCompass />,
      color: '#22c55e',
      description: 'Advanced molecular and protein structure exploration tools',
      content: {
        subtitle: 'Comprehensive Structural Analysis Suite',
        features: [
          'Molecular structure visualization and analysis',
          'Protein sequence and structure exploration',
          'Polymer characterization and analysis',
          'Interactive 3D structure viewers',
          'Chemical property calculations',
          'Structure-activity relationship analysis'
        ],
        tools: [
          { name: 'Molecule Explorer', description: 'SMILES-based molecular analysis', icon: <FaMolecule /> },
          { name: 'Protein Explorer', description: 'UniProt & PDB structure analysis', icon: <FaDna /> },
          { name: 'Polymer Explorer', description: 'Polymer structure characterization', icon: <FaAtom /> }
        ],
        infographics: {
          analysisMetrics: [
            { value: '2M+', label: 'Structures', color: '#22c55e' },
            { value: '50ms', label: 'Avg Time', color: '#16a34a' },
            { value: '15+', label: 'Properties', color: '#15803d' }
          ],
          capabilityFlow: [
            {
              title: 'Input Structure',
              description: 'SMILES, PDB, or draw molecular structures',
              icon: <FaMolecule />,
              color: '#22c55e'
            },
            {
              title: 'Analyze Properties',
              description: 'Calculate physicochemical properties',
              icon: <FaAtom />,
              color: '#16a34a'
            },
            {
              title: 'Visualize Results',
              description: '3D visualization and interactive exploration',
              icon: <FaEye />,
              color: '#15803d'
            }
          ]
        }
      }
    },
    {
      id: 'identification',
      title: 'Molecular Identification',
      icon: <FaBalanceScale />,
      color: '#22c55e',
      description: 'Powerful identification and similarity search capabilities',
      content: {
        subtitle: 'Advanced Molecular Recognition',
        features: [
          'Similarity search algorithms',
          'Substructure matching',
          'Chemical fingerprint comparison',
          'Cross-database identification',
          'Tanimoto coefficient calculations',
          'Multi-parameter similarity analysis'
        ],
        databases: [
          { name: 'ChEMBL', compounds: '2M+', icon: <FaDatabase /> },
          { name: 'PubChem', compounds: '100M+', icon: <FaDatabase /> },
          { name: 'Custom Sets', compounds: 'Unlimited', icon: <FaUsers /> }
        ],
        infographics: {
          searchMetrics: [
            { label: 'Search Speed', percentage: 95, color: '#ff9a56', icon: <FaRocket /> },
            { label: 'Match Accuracy', percentage: 92, color: '#ffad7a', icon: <FaFingerprint /> },
            { label: 'Database Coverage', percentage: 89, color: '#ffc09f', icon: <FaGlobeAmericas /> }
          ],
          databaseStats: [
            { value: '100M+', label: 'Compounds', color: '#ff9a56' },
            { value: '3', label: 'Databases', color: '#ffad7a' },
            { value: '<1s', label: 'Search Time', color: '#ffc09f' }
          ]
        }
      }
    },
    {
      id: 'optimization',
      title: 'Drug Optimization',
      icon: <FaFlask />,
      color: '#22c55e',
      description: 'AI-powered molecular optimization and generation',
      content: {
        subtitle: 'Next-Generation Drug Design',
        features: [
          'ADMET property prediction',
          'BRICS-based molecular generation',
          'LSTM neural network optimization',
          'Structure-based drug design',
          'Lead compound optimization',
          'Synthetic accessibility prediction'
        ],
        algorithms: [
          { name: 'BRICS', description: 'Fragment-based generation', accuracy: '85%' },
          { name: 'LSTM', description: 'Neural network optimization', accuracy: '92%' },
          { name: 'ADMET-AI', description: 'Property prediction', accuracy: '88%' }
        ],
        infographics: {
          optimizationMetrics: [
            { label: 'ADMET Accuracy', percentage: 88, color: '#a8edea', icon: <FaShieldAlt /> },
            { label: 'Generation Speed', percentage: 92, color: '#d299c2', icon: <FaMagic /> },
            { label: 'Drug-likeness', percentage: 85, color: '#fed6e3', icon: <FaBullseye /> }
          ],
          mlProgress: [
            {
              title: 'Input Molecule',
              description: 'Provide target molecule or constraints',
              icon: <FaMolecule />,
              color: '#a8edea'
            },
            {
              title: 'AI Processing',
              description: 'LSTM and BRICS algorithms optimize',
              icon: <FaBrain />,
              color: '#d299c2'
            },
            {
              title: 'Optimized Output',
              description: 'Enhanced molecules with better properties',
              icon: <FaArrowUp />,
              color: '#fed6e3'
            }
          ]
        }
      }
    },
    {
      id: 'activities',
      title: 'Workflow Management',
      icon: <FaTasks />,
      color: '#fed6e3',
      description: 'Comprehensive project and workflow management',
      content: {
        subtitle: 'Streamlined Research Management',
        features: [
          'Project workflow tracking',
          'Task dependency management',
          'Progress visualization',
          'Result organization',
          'Collaborative workspaces',
          'Automated reporting'
        ],
        capabilities: [
          'Multi-step workflow design',
          'Real-time progress tracking',
          'Automated quality checks',
          'Export and sharing tools'
        ],
        infographics: {
          workflowMetrics: [
            { label: 'Project Efficiency', percentage: 96, color: '#fed6e3', icon: <FaTasks /> },
            { label: 'Collaboration', percentage: 89, color: '#d299c2', icon: <FaUsers /> },
            { label: 'Automation', percentage: 94, color: '#a8edea', icon: <FaSyncAlt /> }
          ],
          projectStats: [
            { value: '1000+', label: 'Projects', color: '#fed6e3' },
            { value: '24/7', label: 'Monitoring', color: '#d299c2' },
            { value: '99%', label: 'Success Rate', color: '#a8edea' }
          ]
        }
      }
    }
  ];

  const quickActions = [
    {
      title: 'Start Chat with Master Agent',
      description: 'Begin a conversation with our AI assistant',
      icon: <FaUserSecret />,
      color: '#22c55e',
      route: '/chatbot'
    },
    {
      title: 'Explore Molecules',
      description: 'Analyze molecular structures and properties',
      icon: <FaMolecule />,
      color: '#16a34a',
      route: '/explorer'
    },
    {
      title: 'Search Similar Compounds',
      description: 'Find similar molecules in databases',
      icon: <FaSearch />,
      color: '#15803d',
      route: '/identification'
    },
    {
      title: 'Optimize Drug Properties',
      description: 'Use AI to improve molecular properties',
      icon: <FaLab />,
      color: '#14532d',
      route: '/optimization'
    }
  ];

  return (
    <div className="home-container">
      {/* Timeline Components */}
      <div className="timeline-line"></div>
      <div
        className="timeline-progress"
        style={{ height: `${scrollProgress}%` }}
      ></div>

      <div className="timeline-container" ref={timelineRef}>
        {/* Hero Section */}
        <motion.div
          ref={heroRef}
          initial="hidden"
          animate="visible"
          variants={fadeInUpVariantStatic}
          className="hero-section"
        >
          <div className={`timeline-marker ${activeSection === 'hero' ? 'active' : ''}`}></div>
          <div className="hero-content">
            {/* <video autoPlay loop muted className="hero-video">
              <source src="/videos/hero.mp4" type="video/mp4" />
            </video> */}
            <GlassyContainer className="hero-glassy-container">
            <motion.div
              initial={{ opacity: 0, y: -20 }}
              animate={{ opacity: 1, y: 0 }}
              transition={{ duration: 0.8, delay: 0.2 }}
              className="hero-logo"
            >
              <img src="/images/logo_rendered_main.png" alt="Valency Logo" className="main-logo" />
            </motion.div>

            <motion.h1
              initial={{ opacity: 0, y: -20 }}
              animate={{ opacity: 1, y: 0 }}
              transition={{ duration: 0.8, delay: 0.4 }}
              className="hero-title"
            >
              Welcome to Valency
            </motion.h1>

            <motion.p
              initial={{ opacity: 0, y: -20 }}
              animate={{ opacity: 1, y: 0 }}
              transition={{ duration: 0.8, delay: 0.6 }}
              className="hero-subtitle"
            >
              Next-Generation AI-Powered Drug Discovery Platform
            </motion.p>

            <motion.div
              initial={{ opacity: 0, y: 20 }}
              animate={{ opacity: 1, y: 0 }}
              transition={{ duration: 0.8, delay: 0.8 }}
              className="hero-description"
            >
              Harness the power of artificial intelligence to accelerate drug discovery,
              analyze molecular structures, and optimize therapeutic compounds with unprecedented precision.
            </motion.div>
            </GlassyContainer>
          </div>
        </motion.div>

        {/* Quick Actions */}
        <motion.div
          ref={quickActionsRef}
          initial="hidden"
          animate="visible"
          variants={fadeInUpVariantStatic}
          className="quick-actions-section"
        >
          <div className={`timeline-marker ${activeSection === 'quick-actions' ? 'active' : ''}`}></div>
          <h2 className="section-title">Quick Start</h2>
          <div className="quick-actions-grid">
            {quickActions.map((action, index) => (
              <motion.div
                key={index}
                initial="hidden"
                animate="visible"
                variants={fadeInDownVariants}
                custom={index}
                className="quick-action-card"
                onClick={() => navigate(action.route)}
              >
                <GlassyContainer className="action-card-content">
                  <div
                    className="action-icon"
                    style={{ color: action.color }}
                  >
                    {action.icon}
                  </div>
                  <h3 className="action-title">{action.title}</h3>
                  <p className="action-description">{action.description}</p>
                  <div className="action-button">
                    <FaPlay className="play-icon" />
                    Get Started
                  </div>
                </GlassyContainer>
              </motion.div>
            ))}
          </div>
        </motion.div>

        {/* Platform Sections */}
        <div ref={platformRef} className="platform-sections">
          <div className={`timeline-marker ${activeSection === 'platform' ? 'active' : ''}`}></div>
          <h2 className="section-title">Platform Capabilities</h2>

          {platformSections.map((section, index) => (
            <motion.div
              key={section.id}
              initial="hidden"
              animate="visible"
              variants={fadeInLeftVariants}
              custom={index}
              className="platform-section"
            >
              <GlassyContainer className="section-container">
                <div
                  className="section-header"
                  onClick={() => toggleSection(section.id)}
                >
                  <div className="section-header-left">
                    <div
                      className="section-icon"
                      style={{ color: section.color }}
                    >
                      {section.icon}
                    </div>
                    <div className="section-header-text">
                      <h3 className="section-header-title">{section.title}</h3>
                      <p className="section-header-description">{section.description}</p>
                    </div>
                  </div>
                  <div className="section-toggle">
                    {expandedSections.has(section.id) ? <FaChevronUp /> : <FaChevronDown />}
                  </div>
                </div>

                {expandedSections.has(section.id) && (
                  <motion.div
                    initial={{ opacity: 0, height: 0 }}
                    animate={{ opacity: 1, height: 'auto' }}
                    exit={{ opacity: 0, height: 0 }}
                    transition={{ duration: 0.3 }}
                    className="section-content"
                  >
                    <h4 className="content-subtitle">{section.content.subtitle}</h4>

                    <div className="features-grid">
                      <div className="features-list">
                        <h5>Key Features</h5>
                        <ul>
                          {section.content.features.map((feature, idx) => (
                            <li key={idx}>{feature}</li>
                          ))}
                        </ul>
                      </div>

                      {section.content.stats && (
                        <div className="stats-grid">
                          {section.content.stats.map((stat, idx) => (
                            <div key={idx} className="stat-card">
                              <div className="stat-icon">{stat.icon}</div>
                              <div className="stat-value">{stat.value}</div>
                              <div className="stat-label">{stat.label}</div>
                            </div>
                          ))}
                        </div>
                      )}

                      {section.content.agents && (
                        <div className="agents-grid">
                          <h5>Specialized AI Agents</h5>
                          {section.content.agents.map((agent, idx) => (
                            <div key={idx} className="agent-card">
                              <div className="agent-icon">{agent.icon}</div>
                              <div className="agent-info">
                                <div className="agent-name">{agent.name}</div>
                                <div className="agent-specialty">{agent.specialty}</div>
                              </div>
                            </div>
                          ))}
                        </div>
                      )}

                      {section.content.tools && (
                        <div className="tools-grid">
                          <h5>Available Tools</h5>
                          {section.content.tools.map((tool, idx) => (
                            <div key={idx} className="tool-card">
                              <div className="tool-icon">{tool.icon}</div>
                              <div className="tool-info">
                                <div className="tool-name">{tool.name}</div>
                                <div className="tool-description">{tool.description}</div>
                              </div>
                            </div>
                          ))}
                        </div>
                      )}

                      {section.content.databases && (
                        <div className="databases-grid">
                          <h5>Integrated Databases</h5>
                          {section.content.databases.map((db, idx) => (
                            <div key={idx} className="database-card">
                              <div className="database-icon">{db.icon}</div>
                              <div className="database-info">
                                <div className="database-name">{db.name}</div>
                                <div className="database-compounds">{db.compounds} compounds</div>
                              </div>
                            </div>
                          ))}
                        </div>
                      )}

                      {section.content.algorithms && (
                        <div className="algorithms-grid">
                          <h5>AI Algorithms</h5>
                          {section.content.algorithms.map((algo, idx) => (
                            <div key={idx} className="algorithm-card">
                              <div className="algorithm-info">
                                <div className="algorithm-name">{algo.name}</div>
                                <div className="algorithm-description">{algo.description}</div>
                              </div>
                              <div className="algorithm-accuracy">
                                <FaChartLine />
                                {algo.accuracy}
                              </div>
                            </div>
                          ))}
                        </div>
                      )}

                      {section.content.capabilities && (
                        <div className="capabilities-list">
                          <h5>Capabilities</h5>
                          <ul>
                            {section.content.capabilities.map((capability, idx) => (
                              <li key={idx}>{capability}</li>
                            ))}
                          </ul>
                        </div>
                      )}
                    </div>

                    {/* Infographics Section */}
                    <div className="infographics-section">
                      {section.content.infographics && (
                        <>
                          {section.content.infographics.circularStats && (
                            <div className="circular-stats-container">
                              {section.content.infographics.circularStats.map((stat, idx) => (
                                <StatsCircle
                                  key={idx}
                                  value={stat.value}
                                  label={stat.label}
                                  color={stat.color}
                                  maxValue={100}
                                />
                              ))}
                            </div>
                          )}

                          {section.content.infographics.processSteps && (
                            <div className="process-steps-container">
                              <h5>Process Overview</h5>
                              <ProcessFlow steps={section.content.infographics.processSteps} />
                            </div>
                          )}

                          {section.content.infographics.networkGraph && (
                            <div className="network-graph-container">
                              <h5>Agent Network</h5>
                              <NetworkGraph />
                            </div>
                          )}

                          {section.content.infographics.analysisMetrics && (
                            <div className="analysis-metrics-container">
                              <h5>Analysis Metrics</h5>
                              <div className="metrics-grid">
                                {section.content.infographics.analysisMetrics.map((metric, idx) => (
                                  <div key={idx} className="metric-card">
                                    <div className="metric-value" style={{ color: metric.color }}>
                                      {metric.value}
                                    </div>
                                    <div className="metric-label">
                                      {metric.label}
                                    </div>
                                  </div>
                                ))}
                              </div>
                            </div>
                          )}

                          {section.content.infographics.capabilityFlow && (
                            <div className="capability-flow-container">
                              <h5>Capability Flow</h5>
                              <ProcessFlow steps={section.content.infographics.capabilityFlow} />
                            </div>
                          )}

                          {section.content.infographics.searchMetrics && (
                            <div className="search-metrics-container">
                              <h5>Search Performance</h5>
                              <div className="metrics-grid">
                                {section.content.infographics.searchMetrics.map((metric, idx) => (
                                  <div key={idx} className="metric-card">
                                    <div className="metric-icon" style={{ color: metric.color }}>
                                      {metric.icon}
                                    </div>
                                    <div className="metric-info">
                                      <div className="metric-label">{metric.label}</div>
                                      <div className="metric-percentage">{metric.percentage}%</div>
                                    </div>
                                  </div>
                                ))}
                              </div>
                            </div>
                          )}

                          {section.content.infographics.databaseStats && (
                            <div className="database-stats-container">
                              <h5>Database Statistics</h5>
                              <div className="stats-grid">
                                {section.content.infographics.databaseStats.map((stat, idx) => (
                                  <div key={idx} className="stat-card">
                                    <div className="stat-value" style={{ color: stat.color }}>
                                      {stat.value}
                                    </div>
                                    <div className="stat-label">
                                      {stat.label}
                                    </div>
                                  </div>
                                ))}
                              </div>
                            </div>
                          )}

                          {section.content.infographics.optimizationMetrics && (
                            <div className="optimization-metrics-container">
                              <h5>Optimization Metrics</h5>
                              <div className="metrics-grid">
                                {section.content.infographics.optimizationMetrics.map((metric, idx) => (
                                  <div key={idx} className="metric-card">
                                    <div className="metric-icon" style={{ color: metric.color }}>
                                      {metric.icon}
                                    </div>
                                    <div className="metric-info">
                                      <div className="metric-label">{metric.label}</div>
                                      <div className="metric-percentage">{metric.percentage}%</div>
                                    </div>
                                  </div>
                                ))}
                              </div>
                            </div>
                          )}

                          {section.content.infographics.mlProgress && (
                            <div className="ml-progress-container">
                              <h5>Machine Learning Progress</h5>
                              <ProcessFlow steps={section.content.infographics.mlProgress} />
                            </div>
                          )}

                          {section.content.infographics.workflowMetrics && (
                            <div className="workflow-metrics-container">
                              <h5>Workflow Metrics</h5>
                              <div className="metrics-grid">
                                {section.content.infographics.workflowMetrics.map((metric, idx) => (
                                  <div key={idx} className="metric-card">
                                    <div className="metric-icon" style={{ color: metric.color }}>
                                      {metric.icon}
                                    </div>
                                    <div className="metric-info">
                                      <div className="metric-label">{metric.label}</div>
                                      <div className="metric-percentage">{metric.percentage}%</div>
                                    </div>
                                  </div>
                                ))}
                              </div>
                            </div>
                          )}

                          {section.content.infographics.projectStats && (
                            <div className="project-stats-container">
                              <h5>Project Statistics</h5>
                              <div className="stats-grid">
                                {section.content.infographics.projectStats.map((stat, idx) => (
                                  <div key={idx} className="stat-card">
                                    <div className="stat-value" style={{ color: stat.color }}>
                                      {stat.value}
                                    </div>
                                    <div className="stat-label">
                                      {stat.label}
                                    </div>
                                  </div>
                                ))}
                              </div>
                            </div>
                          )}
                        </>
                      )}
                    </div>
                  </motion.div>
                )}
              </GlassyContainer>
            </motion.div>
          ))}
        </div>

        {/* Call to Action */}
        <motion.div
          ref={ctaRef}
          initial="hidden"
          animate="visible"
          variants={fadeInUpVariantStatic}
          className="cta-section"
        >
          <div className={`timeline-marker ${activeSection === 'cta' ? 'active' : ''}`}></div>
          <GlassyContainer className="cta-content">
            <h2 className="cta-title">Ready to Accelerate Your Research?</h2>
            <p className="cta-description">
              Join thousands of researchers using Valency to transform drug discovery
            </p>
            <button
              className="cta-button"
              onClick={() => navigate('/chatbot')}
            >
              <FaUserSecret />
              Start with Master Agent
            </button>
          </GlassyContainer>
        </motion.div>
      </div>
    </div>
  );
};

export default Home;