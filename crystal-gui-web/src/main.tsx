// StrictMode is disabled because it causes double mounting in development mode,
// which breaks MolStar's createPluginUI (cannot call createRoot twice on same container)
import { createRoot } from 'react-dom/client'
import { Provider } from 'react-redux'
import { store } from './store'
import './index.css'
import App from './App'

createRoot(document.getElementById('root')!).render(
  <Provider store={store}>
    <App />
  </Provider>,
)
