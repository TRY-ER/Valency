import React, { useState, useEffect, useRef } from 'react';
import { color, motion } from 'framer-motion';
import { useNavigate } from 'react-router-dom';
import {
  FaCompass,
  FaDesktop,
  FaRoute,
  FaBalanceScale,
  FaThumbsUp,
  FaFlask,
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
  FaArrowUp,
  FaChevronLeft,
  FaChevronRight,
  FaClock
} from 'react-icons/fa';
import GlassyContainer from '../components/glassy_container/gc';
import {
  fadeInUpVariantStatic,
  fadeInDownVariants,
  fadeInLeftVariants,
  fadeInRightVariants,
  scaleInVariants,
  slideInVariants,
  bounceInVariants,
  rotateInVariants,
  containerVariants,
  itemVariants,
  hoverVariants,
  cardVariants,
  listItemVariants,
  expandableVariants,
  fadeInUpVariants
} from '../components/animations/framerAnim';
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
          <div className="process-step-left">
            <div className="process-step-number-container">
              <div className="process-step-number">{index + 1}</div>
              {index < steps.length - 1 && (
                <div className="process-step-line"></div>
              )}
            </div>
          </div>
          <div className="process-step-right">
            <div className="process-step-header">
              <div className="process-step-icon" style={{ color: step.color }}>
                {step.icon}
              </div>
              <h6 className="process-step-title">{step.title}</h6>
            </div>
            <div className="process-step-content">
              <p className="process-step-description">{step.description}</p>
            </div>
          </div>
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
  const quickActionsGridRef = useRef(null); // Added ref for the grid

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

  const scrollQuickActions = (direction) => {
    if (quickActionsGridRef.current) {
      const scrollAmount = direction === 'left' ? -300 : 300; // Adjust scroll amount as needed
      quickActionsGridRef.current.scrollBy({ left: scrollAmount, behavior: 'smooth' });
    }
  };

  // Enhanced platform sections with infographics
  const platformSections = [
    {
      id: 'overview',
      title: 'Platform Overview',
      icon: <FaLightbulb />,
      color: '#22c55e',
      description: 'Agentic drug discovery and molecular analysis platform',
      content: {
        subtitle: 'Advanced AI-Powered Drug Discovery Toolkit',
        features: [
          'Scientific Tool Based AI Agent Collaborative System',
          'Integrated molecular analysis and visualization',
          'AI-driven optimization and generation algorithms',
          'Comprehensive database integration (ChEMBL, PubChem, UniProt)',
          'Real-time ADMET property prediction',
          'Protein structure analysis and prediction',
          'Interactive 3D molecular visualization'
        ],
        stats: [
          { label: 'AI Agents', value: '10+', icon: <FaRobot />, color: '#0edf5b' },
          { label: 'Scientific Databases', value: '5+', icon: <FaDatabase />, color: '#0fa346' },
          { label: 'Research Ideas', value: '∞', icon: <FaGraduationCap />, color: '#04722d' },
        ],
        infographics: {
          platformMetrics: [
            { label: 'System Reliability', percentage: 99, color: 'var(--color-success)', icon: <FaShieldAlt /> },
            { label: 'AI Accuracy', percentage: 94, color: 'var(--color-accent)', icon: <FaBrain /> },
            { label: 'User Satisfaction', percentage: 97, color: 'var(--color-success)', icon: <FaThumbsUp /> }
          ],
          automationWorkflow: [
            {
              title: 'Natural Language Query',
              description: 'Describe your research goal in plain English',
              icon: <FaRobot />,
              color: 'var(--color-accent)'
            },
            {
              title: 'AI Agent Orchestration',
              description: 'Master agent delegates tasks to specialized agents',
              icon: <FaBrain />,
              color: 'var(--color-accent)'
            },
            {
              title: 'Automated Processing',
              description: 'AI agents execute complex workflows autonomously',
              icon: <FaMagic />,
              color: 'var(--color-accent)'
            },
            {
              title: 'Tool Monitoring',
              description: 'Real-time tracking of analysis and results from various scientific algorithms',
              icon: <FaDesktop />,
              color: 'var(--color-success)'
            },
            {
              title: 'Intelligent Results',
              description: 'Comprehensive insights delivered automatically',
              icon: <FaChartLine />,
              color: 'var(--color-success)'
            }

          ],
          manualWorkflow: [
            {
              title: 'Data Input',
              description: 'Input correspoding scientific data of protines, DNAs, Drugs, etc.',
              icon: <FaUpload />,
              color: 'var(--color-success)'
            },
            {
              title: 'Tool Selection',
              description: 'Choose specific analysis tools and parameters',
              icon: <FaCog />,
              color: 'var(--color-success)'
            },
            {
              title: 'Manual Analysis',
              description: 'Execute step-by-step analysis with full control',
              icon: <FaPlay />,
              color: 'var(--color-success)'
            },
            {
              title: 'Relevant Redirections',
              description: 'Redirect to relevant tools and resources with ease following website structure.',
              icon: <FaRoute />,
              color: 'var(--color-success)'
            },
            {
              title: 'Custom Results',
              description: 'Tailored visualizations and detailed reports of candidates',
              icon: <FaChartPie />,
              color: 'var(--color-success)'
            }
          ]
        }
      }
    },
    {
      id: 'master-agent',
      title: 'Master Agent System',
      icon: <FaRobot />,
      color: '#22c55e',
      description: 'Intelligent agentic orchestrator for complex drug discovery workflows',
      content: {
        subtitle: 'Your AI Research Assistant',
        features: [
          'Natural language query processing',
          'Multi-agent coordination and delegation',
          'Context-aware task routing',
          'Specialized sub-agents for different domains',
          'Completely functional tool interaction within chat interface',
          'Real-time workflow management',
          'Interactive conversational interface'
        ],
        agents: [
          { name: 'Protein Agent', specialty: 'Protein analysis & structure', icon: <FaDna /> },
          { name: 'Drug Agent', specialty: 'Database searches & exploration', icon: <FaSearch /> },
          { name: 'Drug Optimization Agent', specialty: 'ADMET prediction & molecular generation', icon: <FaFlask /> }
        ]
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
          {
            name: 'Molecule Explorer',
            description: 'SMILES-based molecular analysis',
            icon: <FaMolecule />,
            route: '/explorer',
            buttonText: 'Explore Molecules'
          },
          {
            name: 'Protein Explorer',
            description: 'UniProt & PDB structure analysis',
            icon: <FaDna />,
            route: '/explorer/proe',
            buttonText: 'Analyze Proteins'
          },
          {
            name: 'Polymer Explorer',
            description: 'Polymer structure Analysis',
            icon: <FaAtom />,
            route: 'explorer/polye',
            buttonText: 'Study Polymers'
          }
        ],
        infographics: {
          analysisMetrics: [
            { value: '100M+', label: 'Structures', color: '#24e76b', icon: <FaDatabase /> },
            { value: '50ms', label: 'Avg Time', color: '#05ad42', icon: <FaClock /> },
            { value: '50+', label: 'Properties', color: '#016b28', icon: <FaAtom /> }
          ],
          capabilityFlow: [
            {
              title: 'Input Structure',
              description: 'PDB, UniProt, SMILES, or PSMILES input',
              icon: <FaMolecule />,
              color: '#24e76b'
            },
            {
              title: 'Analyze Properties',
              description: 'Calculate physicochemical and ADMET properties',
              icon: <FaAtom />,
              color: '#05ad42'
            },
            {
              title: 'Visualize Results',
              description: '2D Chemical and 3D Atomic interactive exploration',
              icon: <FaEye />,
              color: '#016b28'
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
          'Multi-parameter similarity analysis',
          'Target and activity analysis',
        ],
        databases: [
          { name: 'ChEMBL', compounds: '1.6M+ Compounds', icon: <FaDatabase />, route: '/identification/chembl', buttonText: 'Search ChEMBL' },
          { name: 'PubChem', compounds: '100M+ Compounds', icon: <FaDatabase />, route: '/identification/pubchem', buttonText: 'Search PubChem' },
          { name: 'UniProt', compounds: '150M+ Candidates', icon: <FaDatabase />, route: '/explorer/uniprot', buttonText: 'Explore UniProt' },
        ],
        infographics: {
          databaseStats: [
            { value: '200M+', label: 'Compounds and Candidates', color: "#03d450", icon: <FaDatabase /> },
            { value: '3+', label: 'Databases', color: '#009938', icon: <FaDatabase /> },
            { value: '<1s', label: 'Search Time', color: '#01501e', icon: <FaClock /> }
          ]
        }
      }
    },
    {
      id: 'optimization',
      title: 'Drug Optimization',
      icon: <FaFlask />,
      color: '#22c55e',
      description: 'AI-powered molecular generation and optimization',
      content: {
        subtitle: 'Next-Generation Drug Design',
        features: [
          'ADMET property prediction',
          'BRICS-based molecular generation',
          'LSTM molecular candidates generation',
          'Structure-based drug design',
          'Lead compound optimization',
        ],
        algorithms: [
          { name: 'BRICS', description: 'Fragment-based generation', icon: <FaMolecule />, route: 'optimization', buttonText: 'Use BRICS' },
          { name: 'LSTM', description: 'Neural network optimization', icon: <FaBrain />, route: 'optimization/lstm', buttonText: 'Use LSTM' },
          { name: 'ADMET-AI', description: 'Property prediction', icon: <FaShieldAlt />, route: 'optimization/admet', buttonText: 'Use ADMET-AI' }
        ],
        infographics: {
          mlProgress: [
            {
              title: 'Input Molecule',
              description: 'Provide target molecule or constraints',
              icon: <FaMolecule />,
              color: '#22c55e'
            },
            {
              title: 'AI Processing',
              description: 'LSTM and BRICS algorithms for generation tasks',
              icon: <FaBrain />,
              color: '#16a34a'
            },
            {
              title: 'Optimized Output',
              description: 'Filter candidates based on ADMET properties',
              icon: <FaArrowUp />,
              color: '#15803d'
            }
          ]
        }
      }
    },
    {
      id: 'activities',
      title: 'Activity Management',
      icon: <FaTasks />,
      color: '#22c55e',
      description: 'Comprehensive activity and output management',
      content: {
        subtitle: 'Streamlined Research with Agent Interaction Management',
        features: [
          'Activity Tracking',
          'Agent and Tool Call Statistics',
          'User Specific Activity Logs',
        ],
        capabilities: [
          { name: 'Multi-step workflow design', description: 'Create complex research workflows maintained within activity.', icon: <FaCodeBranch /> },
          { name: 'Fail-safe Activity Log', description: 'Even after closing the browser the activity is stored with it\'s progress.', icon: <FaCheckCircle /> },
          { name: 'Export and using tools', description: 'Activity shows the tool usage as counts to ensure conditionally show tools and preserver resources.', icon: <FaUpload /> }
        ]
      }
    }
  ];

  const quickActions = [
    {
      title: 'Start Chat with Master Agent',
      description: 'Begin a conversation with our Automation Agent',
      icon: <FaRobot className='action-icon-object' />,
      route: '/chatbot'
    },
    {
      title: 'Explore Molecules',
      description: 'Analyze molecular structures and properties',
      icon: <FaMolecule className='action-icon-object' />,
      route: '/explorer'
    },
    {
      title: 'Search Similar Compounds',
      description: 'Find similar molecules in databases',
      icon: <FaSearch className='action-icon-object' />,
      route: '/identification'
    },
    {
      title: 'Optimize Drug Properties',
      description: 'Use AI to generate and optimize compounds',
      icon: <FaLab className='action-icon-object' />,
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
          <div className={`timeline-marker active`}></div>
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
                The Agentic Drug Discovery Platform
              </motion.p>

              <motion.div
                initial={{ opacity: 0, y: 20 }}
                animate={{ opacity: 1, y: 0 }}
                transition={{ duration: 0.8, delay: 0.8 }}
                className="hero-description"
              >
                Harness the power of artificial intelligence to accelerate drug discovery,
                analyze molecular structures, and optimize therapeutic compounds with unprecedented precision just with your prompt.
              </motion.div>
            </GlassyContainer>
          </div>
        </motion.div>

        {/* Quick Actions */}
        <motion.div
          ref={quickActionsRef}
          initial="hidden"
          animate="visible"
          variants={containerVariants}
          className="quick-actions-section"
        >
          <div className={`timeline-marker active`}></div>
          <motion.h2
            className="section-title"
            variants={scaleInVariants}
            initial="hidden"
            animate="visible"
          >
            Quick Start
          </motion.h2>
          <div className="quick-actions-grid-wrapper">
            <motion.button
              className="quick-actions-nav-arrow quick-actions-nav-left"
              onClick={() => scrollQuickActions('left')}
              aria-label="Scroll left"
              whileHover={{ scale: 1.1, x: -5 }}
              whileTap={{ scale: 0.9 }}
              transition={{ duration: 0.2 }}
            >
              <FaChevronLeft />
            </motion.button>
            <motion.div
              className="quick-actions-grid"
              ref={quickActionsGridRef}
              variants={containerVariants}
              initial="hidden"
              animate="visible"
            >
              {quickActions.map((action, index) => (
                <motion.div
                  key={index}
                  variants={cardVariants}
                  custom={index}
                  className="quick-action-card"
                  onClick={() => navigate(action.route)}
                  whileHover="hover"
                  whileTap={{ scale: 0.95 }}
                  transition={{ duration: 0.2 }}
                >
                  <GlassyContainer className="action-card-content">
                    <motion.div
                      className="action-icon"
                      variants={rotateInVariants}
                      custom={index}
                    >
                      {action.icon}
                    </motion.div>
                    <motion.h3
                      className="action-title"
                      variants={itemVariants}
                    >
                      {action.title}
                    </motion.h3>
                    <motion.p
                      className="action-description"
                      variants={itemVariants}
                    >
                      {action.description}
                    </motion.p>
                    <motion.div
                      className="action-button"
                      variants={bounceInVariants}
                      custom={index}
                    >
                      <FaPlay className="play-icon" />
                      Get Started
                    </motion.div>
                  </GlassyContainer>
                </motion.div>
              ))}
            </motion.div>
            <motion.button
              className="quick-actions-nav-arrow quick-actions-nav-right"
              onClick={() => scrollQuickActions('right')}
              aria-label="Scroll right"
              whileHover={{ scale: 1.1, x: 5 }}
              whileTap={{ scale: 0.9 }}
              transition={{ duration: 0.2 }}
            >
              <FaChevronRight />
            </motion.button>
          </div>
        </motion.div>

        {/* Platform Sections */}
        <motion.div
          ref={platformRef}
          className="platform-sections"
          variants={containerVariants}
          initial="hidden"
          animate="visible"
        >
          <div className={`timeline-marker active`}></div>
          <motion.h2
            className="section-title"
            variants={scaleInVariants}
            initial="hidden"
            animate="visible"
          >
            Platform Capabilities
          </motion.h2>

          {platformSections.map((section, index) => (
            <motion.div
              key={section.id}
              variants={slideInVariants}
              custom={index}
              className="platform-section"
            >
              <GlassyContainer className="section-container">
                <motion.div
                  className="section-header"
                  onClick={() => toggleSection(section.id)}
                  whileHover={{ scale: 1.02 }}
                  whileTap={{ scale: 0.98 }}
                  transition={{ duration: 0.2 }}
                >
                  <div className="section-header-left">
                    <motion.div
                      className="section-icon"
                      style={{ color: section.color }}
                      variants={rotateInVariants}
                      initial="hidden"
                      animate="visible"
                      custom={index}
                    >
                      {section.icon}
                    </motion.div>
                    <div className="section-header-text">
                      <motion.h3
                        className="section-header-title"
                        variants={fadeInRightVariants}
                        initial="hidden"
                        animate="visible"
                        custom={index}
                      >
                        {section.title}
                      </motion.h3>
                      <motion.p
                        className="section-header-description"
                        variants={fadeInRightVariants}
                        initial="hidden"
                        animate="visible"
                        custom={index + 0.5}
                      >
                        {section.description}
                      </motion.p>
                    </div>
                  </div>
                  <motion.div
                    className="section-toggle"
                    animate={{
                      rotate: expandedSections.has(section.id) ? 180 : 0
                    }}
                    transition={{ duration: 0.3 }}
                  >
                    <FaChevronDown />
                  </motion.div>
                </motion.div>

                {expandedSections.has(section.id) && (
                  <motion.div
                    variants={expandableVariants}
                    initial="collapsed"
                    animate="expanded"
                    exit="exit"
                    className="section-content"
                  >
                    <motion.h4
                      className="content-subtitle"
                      variants={scaleInVariants}
                      initial="hidden"
                      animate="visible"
                    >
                      {section.content.subtitle}
                    </motion.h4>

                    <motion.div
                      className="features-grid"
                      variants={containerVariants}
                      initial="hidden"
                      animate="visible"
                    >
                      <motion.div
                        className="features-list"
                        variants={itemVariants}
                      >
                        <h5>Key Features</h5>
                        <ul>
                          {section.content.features.map((feature, idx) => (
                            <motion.li
                              key={idx}
                              variants={listItemVariants}
                              custom={idx}
                              initial="hidden"
                              animate="visible"
                            >
                              {feature}
                            </motion.li>
                          ))}
                        </ul>
                      </motion.div>

                      {section.content.stats && (
                        <motion.div
                          className="stats-grid"
                          variants={containerVariants}
                          initial="hidden"
                          animate="visible"
                        >
                          {section.content.stats.map((stat, idx) => (
                            <motion.div
                              key={idx}
                              className="stat-card"
                              variants={cardVariants}
                              custom={idx}
                              whileHover="hover"
                            >
                              <motion.div
                                className="stat-icon"
                                style={{ color: `${stat.color}` }}
                                variants={bounceInVariants}
                                custom={idx}
                              >
                                {stat.icon}
                              </motion.div>
                              <motion.div
                                className="stat-value"
                                variants={scaleInVariants}
                                custom={idx}
                              >
                                {stat.value}
                              </motion.div>
                              <motion.div
                                className="stat-label"
                                variants={fadeInUpVariants}
                                custom={idx}
                              >
                                {stat.label}
                              </motion.div>
                            </motion.div>
                          ))}
                        </motion.div>
                      )}

                      {section.content.agents && (
                        <motion.div
                          className="agents-grid"
                          variants={containerVariants}
                          initial="hidden"
                          animate="visible"
                        >
                          <h5>Specialized AI Agents</h5>
                          {section.content.agents.map((agent, idx) => (
                            <motion.div
                              key={idx}
                              className="agent-card"
                              variants={cardVariants}
                              custom={idx}
                              whileHover="hover"
                            >
                              <motion.div
                                className="agent-icon"
                                variants={rotateInVariants}
                                custom={idx}
                              >
                                {agent.icon}
                              </motion.div>
                              <div className="agent-info">
                                <motion.div
                                  className="agent-name"
                                  variants={fadeInUpVariants}
                                  custom={idx}
                                >
                                  {agent.name}
                                </motion.div>
                                <motion.div
                                  className="agent-specialty"
                                  variants={fadeInUpVariants}
                                  custom={idx + 0.2}
                                >
                                  {agent.specialty}
                                </motion.div>
                              </div>
                            </motion.div>
                          ))}
                          
                          {/* Add chatbot button for master-agent section only */}
                          {section.id === 'master-agent' && (
                            <motion.div
                              className="master-agent-button-container"
                              variants={containerVariants}
                              initial="hidden"
                              animate="visible"
                              style={{ marginTop: '25px', textAlign: 'center' }}
                            >
                              <motion.button
                                className="master-agent-redirect-button"
                                onClick={() => navigate('/chatbot')}
                                variants={scaleInVariants}
                                whileHover={{
                                  backgroundColor: '#16a34a',
                                  y: -2,
                                  boxShadow: '0 6px 12px rgba(0,0,0,0.15)',
                                  scale: 1.05
                                }}
                                whileTap={{ scale: 0.95 }}
                                style={{
                                  padding: '12px 24px',
                                  backgroundColor: '#22c55e',
                                  color: 'white',
                                  border: 'none',
                                  borderRadius: '8px',
                                  cursor: 'pointer',
                                  fontSize: '16px',
                                  fontWeight: '600',
                                  transition: 'all 0.2s ease',
                                  boxShadow: '0 4px 8px rgba(0,0,0,0.1)',
                                  display: 'inline-flex',
                                  alignItems: 'center',
                                  gap: '8px'
                                }}
                              >
                                <FaRobot />
                                Start Chatting with Master Agent
                              </motion.button>
                            </motion.div>
                          )}
                        </motion.div>
                      )}

                      {section.content.tools && (
                        <motion.div
                          className="tools-grid"
                          variants={containerVariants}
                          initial="hidden"
                          animate="visible"
                        >
                          <h5>Available Tools</h5>
                          {section.content.tools.map((tool, idx) => (
                            <motion.div
                              key={idx}
                              className="tool-card"
                              style={{ display: 'flex', alignItems: 'center', justifyContent: 'space-between' }}
                              variants={cardVariants}
                              custom={idx}
                              whileHover="hover"
                            >
                              <div style={{ display: 'flex', alignItems: 'center', flex: 1 }}>
                                <motion.div
                                  className="tool-icon"
                                  style={{ marginRight: '16px' }}
                                  variants={bounceInVariants}
                                  custom={idx}
                                >
                                  {tool.icon}
                                </motion.div>
                                <div className="tool-info">
                                  <motion.div
                                    className="tool-name"
                                    variants={fadeInRightVariants}
                                    custom={idx}
                                  >
                                    {tool.name}
                                  </motion.div>
                                  <motion.div
                                    className="tool-description"
                                    variants={fadeInRightVariants}
                                    custom={idx + 0.2}
                                  >
                                    {tool.description}
                                  </motion.div>
                                </div>
                              </div>
                              {tool.route && tool.buttonText && (
                                <motion.button
                                  className="tool-redirect-button"
                                  onClick={() => navigate(tool.route)}
                                  style={{
                                    marginLeft: '20px',
                                    padding: '8px 16px',
                                    backgroundColor: '#22c55e',
                                    color: 'white',
                                    border: 'none',
                                    borderRadius: '6px',
                                    cursor: 'pointer',
                                    fontSize: '14px',
                                    fontWeight: '600',
                                    transition: 'all 0.2s ease',
                                    boxShadow: '0 2px 4px rgba(0,0,0,0.1)',
                                    flexShrink: 0
                                  }}
                                  whileHover={{
                                    backgroundColor: '#16a34a',
                                    y: -1,
                                    boxShadow: '0 4px 8px rgba(0,0,0,0.15)',
                                    scale: 1.05
                                  }}
                                  whileTap={{ scale: 0.95 }}
                                  variants={scaleInVariants}
                                  custom={idx}
                                >
                                  {tool.buttonText}
                                </motion.button>
                              )}
                            </motion.div>
                          ))}
                        </motion.div>
                      )}

                      {section.content.databases && (
                        <div className="databases-grid">
                          <h5>Integrated Databases</h5>
                          {section.content.databases.map((db, idx) => (
                            <div key={idx} className="database-card" style={{ display: 'flex', alignItems: 'center', justifyContent: 'space-between' }}>
                              <div style={{ display: 'flex', alignItems: 'center', flex: 1 }}>
                                <div className="database-icon" style={{ marginRight: '16px' }}>{db.icon}</div>
                                <div className="database-info">
                                  <div className="database-name">{db.name}</div>
                                  <div className="database-compounds">{db.compounds}</div>
                                </div>
                              </div>
                              {db.route && db.buttonText && (
                                <button
                                  className="database-redirect-button"
                                  onClick={() => navigate(db.route)}
                                  style={{
                                    marginLeft: '20px',
                                    padding: '8px 16px',
                                    backgroundColor: '#22c55e',
                                    color: 'white',
                                    border: 'none',
                                    borderRadius: '6px',
                                    cursor: 'pointer',
                                    fontSize: '14px',
                                    fontWeight: '600',
                                    transition: 'all 0.2s ease',
                                    boxShadow: '0 2px 4px rgba(0,0,0,0.1)',
                                    flexShrink: 0
                                  }}
                                  onMouseOver={(e) => {
                                    e.target.style.backgroundColor = '#16a34a';
                                    e.target.style.transform = 'translateY(-1px)';
                                    e.target.style.boxShadow = '0 4px 8px rgba(0,0,0,0.15)';
                                  }}
                                  onMouseOut={(e) => {
                                    e.target.style.backgroundColor = '#22c55e';
                                    e.target.style.transform = 'translateY(0)';
                                    e.target.style.boxShadow = '0 2px 4px rgba(0,0,0,0.1)';
                                  }}
                                >
                                  {db.buttonText}
                                </button>
                              )}
                            </div>
                          ))}
                        </div>
                      )}

                      {section.content.algorithms && (
                        <div className="algorithms-grid">
                          <h5>AI Algorithms</h5>
                          {section.content.algorithms.map((algo, idx) => (
                            <div key={idx} className="algorithm-card" style={{ display: 'flex', alignItems: 'center', justifyContent: 'space-between' }}>
                              <div style={{ display: 'flex', alignItems: 'center', flex: 1 }}>
                                <div className="algorithm-icon" style={{ marginRight: '16px' }}>{algo.icon}</div>
                                <div className="algorithm-info">
                                  <div className="algorithm-name">{algo.name}</div>
                                  <div className="algorithm-description">{algo.description}</div>
                                </div>
                              </div>
                              {algo.route && algo.buttonText && (
                                <button
                                  className="algorithm-redirect-button"
                                  onClick={() => navigate(algo.route)}
                                  style={{
                                    marginLeft: '20px',
                                    padding: '8px 16px',
                                    backgroundColor: '#22c55e',
                                    color: 'white',
                                    border: 'none',
                                    borderRadius: '6px',
                                    cursor: 'pointer',
                                    fontSize: '14px',
                                    fontWeight: '600',
                                    transition: 'all 0.2s ease',
                                    boxShadow: '0 2px 4px rgba(0,0,0,0.1)',
                                    flexShrink: 0
                                  }}
                                  onMouseOver={(e) => {
                                    e.target.style.backgroundColor = '#16a34a';
                                    e.target.style.transform = 'translateY(-1px)';
                                    e.target.style.boxShadow = '0 4px 8px rgba(0,0,0,0.15)';
                                  }}
                                  onMouseOut={(e) => {
                                    e.target.style.backgroundColor = '#22c55e';
                                    e.target.style.transform = 'translateY(0)';
                                    e.target.style.boxShadow = '0 2px 4px rgba(0,0,0,0.1)';
                                  }}
                                >
                                  {algo.buttonText}
                                </button>
                              )}
                            </div>
                          ))}
                        </div>
                      )}

                      {section.content.capabilities && (
                        <div className="capabilities-grid">
                          <h5>Advanced Capabilities</h5>
                          {section.content.capabilities.map((capability, idx) => (
                            <div key={idx} className="capability-card">
                              <div className="capability-icon">{capability.icon}</div>
                              <div className="capability-info">
                                <div className="capability-name">{capability.name}</div>
                                <div className="capability-description">{capability.description}</div>
                              </div>
                            </div>
                          ))}
                        </div>
                      )}

                      {/* Add activities button for activities section only */}
                      {section.id === 'activities' && (
                        <motion.div
                          className="activities-button-container"
                          variants={containerVariants}
                          initial="hidden"
                          animate="visible"
                          style={{ marginTop: '25px', textAlign: 'center' }}
                        >
                          <motion.button
                            className="activities-redirect-button"
                            onClick={() => navigate('/activities')}
                            variants={scaleInVariants}
                            whileHover={{
                              backgroundColor: '#16a34a',
                              y: -2,
                              boxShadow: '0 6px 12px rgba(0,0,0,0.15)',
                              scale: 1.05
                            }}
                            whileTap={{ scale: 0.95 }}
                            style={{
                              padding: '12px 24px',
                              backgroundColor: '#22c55e',
                              color: 'white',
                              border: 'none',
                              borderRadius: '8px',
                              cursor: 'pointer',
                              fontSize: '16px',
                              fontWeight: '600',
                              transition: 'all 0.2s ease',
                              boxShadow: '0 4px 8px rgba(0,0,0,0.1)',
                              display: 'inline-flex',
                              alignItems: 'center',
                              gap: '8px'
                            }}
                          >
                            <FaTasks />
                            Go to Activities Page
                          </motion.button>
                        </motion.div>
                      )}
                    </motion.div>

                    {/* Infographics Section */}
                    <div className="infographics-section">
                      {section.content.infographics && (
                        <>
                          {/* Automation and Manual Workflows */}
                          {(section.content.infographics.automationWorkflow || section.content.infographics.manualWorkflow) && (
                            <div className="dual-workflow-container">
                              <h5>Research Workflows</h5>
                              <div className="workflows-side-by-side">
                                {section.content.infographics.automationWorkflow && (
                                  <div className="workflow-column">
                                    <div className="workflow-header">
                                      <div className="workflow-header-icon">
                                        <FaRobot style={{ color: 'var(--color-accent)' }} />
                                      </div>
                                      <h6>Automation Workflow</h6>
                                      <p>AI-driven research with natural language queries</p>
                                    </div>
                                    <ProcessFlow steps={section.content.infographics.automationWorkflow} />
                                  </div>
                                )}
                                {section.content.infographics.manualWorkflow && (
                                  <div className="workflow-column">
                                    <div className="workflow-header">
                                      <div className="workflow-header-icon">
                                        <FaCog style={{ color: 'var(--color-success)' }} />
                                      </div>
                                      <h6>Manual Workflow</h6>
                                      <p>Precision control with step-by-step analysis</p>
                                    </div>
                                    <ProcessFlow steps={section.content.infographics.manualWorkflow} />
                                  </div>
                                )}
                              </div>
                            </div>
                          )}

                          {/* Legacy Workflow Steps Support */}
                          {section.content.infographics.workflowSteps && !section.content.infographics.automationWorkflow && !section.content.infographics.manualWorkflow && (
                            <div className="workflow-steps-container">
                              <h5>Platform Workflow</h5>
                              <ProcessFlow steps={section.content.infographics.workflowSteps} />
                            </div>
                          )}

                          {/* Analysis Metrics */}
                          {section.content.infographics.analysisMetrics && (
                            <div className="analysis-metrics-container">
                              <h5>Analysis Statistics</h5>
                              <div className="stats-grid">
                                {section.content.infographics.analysisMetrics.map((stat, idx) => (
                                  <div key={idx} className="stat-card">
                                    {stat.icon && (
                                      <div className="stat-icon" style={{ color: stat.color }}>
                                        {stat.icon}
                                      </div>
                                    )}
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

                          {/* Search Metrics */}
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

                          {/* Database Stats */}
                          {section.content.infographics.databaseStats && (
                            <div className="database-stats-container">
                              <h5>Database Statistics</h5>
                              <div className="stats-grid">
                                {section.content.infographics.databaseStats.map((stat, idx) => (
                                  <div key={idx} className="stat-card">
                                    {stat.icon && (
                                      <div className="stat-icon" style={{ color: stat.color }}>
                                        {stat.icon}
                                      </div>
                                    )}
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

                          {/* Optimization Metrics */}
                          {section.content.infographics.optimizationMetrics && (
                            <div className="optimization-metrics-container">
                              <h5>Optimization Performance</h5>
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

                          {/* Capability Flow */}
                          {section.content.infographics.capabilityFlow && (
                            <div className="capability-flow-container">
                              <h5>Capability Flow</h5>
                              <ProcessFlow steps={section.content.infographics.capabilityFlow} />
                            </div>
                          )}

                          {/* ML Progress */}
                          {section.content.infographics.mlProgress && (
                            <div className="ml-progress-container">
                              <h5>Machine Learning Progress</h5>
                              <ProcessFlow steps={section.content.infographics.mlProgress} />
                            </div>
                          )}

                          {/* Legacy handlers for backward compatibility */}
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
      </motion.div>

      {/* Call to Action */}
      <motion.div
        ref={ctaRef}
        initial="hidden"
        animate="visible"
        variants={fadeInUpVariantStatic}
        className="cta-section"
      >
        <div className={`timeline-marker active`}></div>
        <GlassyContainer className="cta-content">
          <h2 className="cta-title">Ready to Accelerate Your Research?</h2>
          <button
            className="cta-button"
            onClick={() => navigate('/chatbot')}
          >
            <FaRobot />
            Start with Master Agent
          </button>
        </GlassyContainer>
      </motion.div>
    </div>
  </div>);
};

export default Home;